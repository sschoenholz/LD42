using System;
using Unity.Collections;
using Unity.Jobs;
using Unity.Entities;
using UnityEngine;
using UnityEngine.Jobs;
using Unity.Transforms;
using Unity.Mathematics;

namespace ParticleSimulator
{


    public class SimulateFrame : JobComponentSystem
    {
        public const int hash_constant = 1997;

        [Unity.Burst.BurstCompileAttribute]
        struct ComputeVelocityVerletHalfStepJob : IJobParallelFor
        {

            [ReadOnly] public ComponentDataArray<Force> forces;
            [ReadOnly] public float dt;
            [ReadOnly] public Vector3 size;

            public ComponentDataArray<Velocity> velocities;
            public ComponentDataArray<Position> positions;

            public void Execute(int i)
            {
                velocities[i] = new Velocity
                {
                    Value = velocities[i].Value + 0.5f * dt * forces[i].Value
                };

                float3 pos = positions[i].Value;
                pos += velocities[i].Value * dt;

                positions[i] = new Position { Value = pos };
            }
        }

        [Unity.Burst.BurstCompileAttribute]
        struct HashParticlesJob : IJobParallelFor
        {
            // NOTE(schsam): We only hash the velocities so that we can reduce
            // over them to efficiently compute PE. If we're not interested
            // in the PE computation we can just hash the positions.
            [ReadOnly] public ComponentDataArray<Position> positions;
            [ReadOnly] public ComponentDataArray<Velocity> velocities;
            [ReadOnly] public ComponentDataArray<InteractionType> interactions;

            public NativeArray<int> hashes;
            public NativeMultiHashMap<int, int>.Concurrent hash_map;
            public NativeMultiHashMap<int, Position>.Concurrent hash_positions;
            public NativeMultiHashMap<int, Velocity>.Concurrent hash_velocities;
            public NativeMultiHashMap<int, InteractionType>.Concurrent hash_interactions;

            public float cell_size;
            public int3 cells_per_side;
            public Vector3 size;


            public void Execute(int i)
            {
                int x_index = (int)Math.Floor(positions[i].Value.x / cell_size) + 10;
                int y_index = (int)Math.Floor(positions[i].Value.y / cell_size) + 50;
                int z_index = (int)Math.Floor(positions[i].Value.z / cell_size) + 10;

                int hash = (
                    z_index * hash_constant * hash_constant
                    + y_index * hash_constant + x_index);

                hashes[i] = hash;
                hash_map.Add(hash, i);
                hash_positions.Add(hash, positions[i]);
                hash_velocities.Add(hash, velocities[i]);
                hash_interactions.Add(hash, interactions[i]);
            }
        }

        [Unity.Burst.BurstCompileAttribute]
        struct InitializeHashSet : IJobParallelFor
        {
            public NativeArray<int> set_filled;
            public NativeArray<int> set_hash;

            public void Execute(int i)
            {
                set_filled[i] = 0;
                set_hash[i] = 0;
            }
        }

        // NOTE(schsam): I feel like we should iterate through the keys directly.
        [Unity.Burst.BurstCompileAttribute]
        struct BuildHashSetJob : IJob
        {
            [ReadOnly] public NativeArray<int> hashes;
            public NativeArray<int> set_hash;
            public NativeArray<int> set_filled;

            public NativeArray<int> unique_hash_count;

            bool InsertInSet(int hash)
            {
                int id = hash % set_filled.Length;

                while (set_filled[id] != 0)
                {
                    if (set_hash[id] == hash)
                        return false;

                    id++;
                    if (id > set_filled.Length)
                        id = 0;
                }

                set_filled[id] = 1;
                set_hash[id] = hash;

                return true;
            }

            public void Execute()
            {
                unique_hash_count[0] = 0;

                for (int i = 0; i < hashes.Length; i++)
                {
                    if (InsertInSet(hashes[i]))
                    {
                        unique_hash_count[0] = unique_hash_count[0] + 1;
                    }
                }
            }
        }

        [Unity.Burst.BurstCompileAttribute]
        struct ComputeDensityJob : IJobParallelFor
        {
            [ReadOnly] public NativeMultiHashMap<int, int> hash_map;
            [ReadOnly] public NativeMultiHashMap<int, Position> cell_positions;
            [ReadOnly] public NativeMultiHashMap<int, InteractionType> cell_interactions;

            [ReadOnly] public NativeArray<int> set_hash;
            [ReadOnly] public NativeArray<int> set_filled;

            [ReadOnly] public float size;
            [ReadOnly] public float size2;
            [ReadOnly] public float isize9;
            [ReadOnly] public float isize6;

            [ReadOnly] public float3 system_size;
            [ReadOnly] public int3 cells_per_side;
            [ReadOnly] public int num_cells;

            public NativeMultiHashMap<int, float>.Concurrent cell_density;
            public NativeMultiHashMap<int, int>.Concurrent density_hash_map;

            public void PolyKernel(Position p_i, Position p_j, ref float kernel)
            {
                const float iscale = 64.0f * Mathf.PI / 315f;

                float3 dr = p_j.Value - p_i.Value;

                float rsq = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;

                if (rsq < size2)
                    kernel += iscale * isize9 * (size2 - rsq) * (size2 - rsq) * (size2 - rsq);
            }

            public void SpikyKernel(Position p_i, Position p_j, ref float kernel)
            {
                const float iscale = Mathf.PI / 15f;

                float3 dr = p_j.Value - p_i.Value;

                float r = Mathf.Sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);

                if (r < size)
                    kernel += iscale * isize6 * (size - r) * (size - r) * (size - r);
            }

            public float SumDensitiesCell(int i, Position p_i, int cell)
            {
                int j;
                Position p_j;
                InteractionType inter_j;

                float density = 0f;

                NativeMultiHashMapIterator<int> pos_it;
                NativeMultiHashMapIterator<int> hash_it;
                NativeMultiHashMapIterator<int> inter_it;

                if (cell_positions.TryGetFirstValue(cell, out p_j, out pos_it) &&
                    hash_map.TryGetFirstValue(cell, out j, out hash_it) &&
                    cell_interactions.TryGetFirstValue(cell, out inter_j, out inter_it))
                {
                    if (i != j)// && inter_j.Value == Interaction.Particle)
                    {
                        PolyKernel(p_i, p_j, ref density);
                    }

                    while (cell_positions.TryGetNextValue(out p_j, ref pos_it) &&
                           hash_map.TryGetNextValue(out j, ref hash_it) &&
                           cell_interactions.TryGetNextValue(out inter_j, ref inter_it))
                    {
                        if (i != j)// && inter_j.Value == Interaction.Particle)
                        {
                            PolyKernel(p_i, p_j, ref density);
                        }
                    }
                }


                return density;
            }

            public float SumDensities(int i, Position p_i, int3 cell)
            {
                float result = 0f;

                for (int dz = -1; dz <= 1; dz++)
                {
                    int cell_index_z = cell.z + dz;

                    for (int dy = -1; dy <= 1; dy++)
                    {
                        int cell_index_y = cell.y + dy;

                        for (int dx = -1; dx <= 1; dx++)
                        {
                            int cell_index_x = cell.x + dx;

                            int current_cell = cell_index_z * hash_constant * hash_constant +
                                cell_index_y * hash_constant +
                                cell_index_x;

                            result += SumDensitiesCell(i, p_i, current_cell);
                        }
                    }
                }

                return result;
            }

            public void Execute(int id)
            {
                if (set_filled[id] == 0)
                    return;

                int cell_hash = set_hash[id];

                int cell_index_x = cell_hash % hash_constant;
                int cell_index_y = cell_hash % (hash_constant * hash_constant);
                cell_index_y /= hash_constant;
                int cell_index_z = cell_hash / (hash_constant * hash_constant);

                int3 cell_index = new int3(cell_index_x, cell_index_y, cell_index_z);

                int i;
                Position p_i;
                InteractionType inter_i;

                NativeMultiHashMapIterator<int> pos_it;
                NativeMultiHashMapIterator<int> hash_it;
                NativeMultiHashMapIterator<int> inter_it;

                if (cell_positions.TryGetFirstValue(cell_hash, out p_i, out pos_it) &&
                    hash_map.TryGetFirstValue(cell_hash, out i, out hash_it) &&
                    cell_interactions.TryGetFirstValue(cell_hash, out inter_i, out inter_it))
                {

                    //if (inter_i.Value == Interaction.Particle)
                    {
                        float result = SumDensities(i, p_i, cell_index);

                        cell_density.Add(cell_hash, result);
                        density_hash_map.Add(cell_hash, i);
                    }

                    while (cell_positions.TryGetNextValue(out p_i, ref pos_it) &&
                           hash_map.TryGetNextValue(out i, ref hash_it) &&
                           cell_interactions.TryGetNextValue(out inter_i, ref inter_it))
                    {
                        //if (inter_i.Value == Interaction.Particle)
                        {
                            float result = SumDensities(i, p_i, cell_index);

                            cell_density.Add(cell_hash, result);
                            density_hash_map.Add(cell_hash, i);
                        }
                    }
                }

            }
        }

        [Unity.Burst.BurstCompileAttribute]
        struct ComputeForcesJob : IJobParallelFor
        {
            [ReadOnly] public NativeMultiHashMap<int, int> hash_map;
            [ReadOnly] public NativeMultiHashMap<int, Position> cell_positions;
            [ReadOnly] public NativeMultiHashMap<int, Velocity> cell_velocities;
            [ReadOnly] public NativeMultiHashMap<int, InteractionType> cell_interactions;

            [ReadOnly] public NativeMultiHashMap<int, float> cell_densities;
            [ReadOnly] public NativeMultiHashMap<int, int> density_hash_map;

            [ReadOnly] public NativeArray<int> set_hash;
            [ReadOnly] public NativeArray<int> set_filled;

            [ReadOnly] public float size;
            [ReadOnly] public float size2;
            [ReadOnly] public float isize6;
            [ReadOnly] public float isize9;
            [ReadOnly] public float viscosity;
            [ReadOnly] public float repulsion;

            [ReadOnly] public float3 system_size;
            [ReadOnly] public int3 cells_per_side;
            [ReadOnly] public int num_cells;

            [ReadOnly] public float gas_constant;
            [ReadOnly] public float ref_density;

            public NativeMultiHashMap<int, Force>.Concurrent cell_forces;
            public NativeMultiHashMap<int, int>.Concurrent force_hash_map;

            public float3 PolyKernelDerivative(Position p_i, Position p_j)
            {
                const float iscale = 64.0f * Mathf.PI / 315f * 3f;

                float3 dr = p_j.Value - p_i.Value;
                float3 force;

                float rsq = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;

                if (rsq < size2)
                    force = iscale * isize9 * (size2 - rsq) * (size2 - rsq) * dr;
                else
                    force = new float3(0f, 0f, 0f);

                return force;
            }

            public float3 SpikyKernelDerivative(Position p_i, Position p_j)
            {
                const float iscale = Mathf.PI / 15f * 3f;

                float3 dr = p_j.Value - p_i.Value;
                float3 force;

                float r = Mathf.Sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);

                if (r < size)
                    force = iscale * isize6 * (size - r) * (size - r) * dr / r;
                else
                    force = new float3(0f, 0f, 0f);

                return force;
            }

            public float ViscosityKernel(Position p_i, Position p_j)
            {
                const float iscale = 45f / Mathf.PI;

                float3 dr = p_j.Value - p_i.Value;
                float r = Mathf.Sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);

                if (r < size)
                    return iscale * isize6 * (size - r);

                return 0f;
            }

            public float3 ObjectKernel(Position p_i, Position p_j)
            {
                float3 dr = p_j.Value - p_i.Value;

                float r = Mathf.Sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);

                const float iscale = Mathf.PI / 15f * 3f;

                if (r < size)
                    return iscale / size / size * (size - r) * dr / r;

                return new float3(0f, 0f, 0f);
            }

            public float GetDensity(int cell, int i)
            {
                NativeMultiHashMapIterator<int> den_it;
                NativeMultiHashMapIterator<int> hash_it;

                float d;
                int p;

                if (cell_densities.TryGetFirstValue(cell, out d, out den_it) &&
                    density_hash_map.TryGetFirstValue(cell, out p, out hash_it))
                {
                    if (p == i)
                        return d;

                    while (cell_densities.TryGetNextValue(out d, ref den_it) &&
                          density_hash_map.TryGetNextValue(out p, ref hash_it))
                    {
                        if (p == i)
                            return d;
                    }
                }

                return -1f;
            }

            public float3 SumForcesCell(int i, Position p_i, Velocity v_i, InteractionType inter_i, float pr_i, int cell)
            {
                int j;
                Position p_j;
                Velocity v_j;
                InteractionType inter_j;

                float3 force = new float3(0f, 0f, 0f);

                NativeMultiHashMapIterator<int> pos_it;
                NativeMultiHashMapIterator<int> vel_it;
                NativeMultiHashMapIterator<int> inter_it;
                NativeMultiHashMapIterator<int> hash_it;

                if (cell_positions.TryGetFirstValue(cell, out p_j, out pos_it) &&
                    hash_map.TryGetFirstValue(cell, out j, out hash_it) &&
                    cell_velocities.TryGetFirstValue(cell, out v_j, out vel_it) &&
                    cell_interactions.TryGetFirstValue(cell, out inter_j, out inter_it))
                {
                    if (i != j)
                    {
                        float d_j = GetDensity(cell, j);
                        float pr_j = gas_constant * (d_j - ref_density);

                        if (inter_i.Value == Interaction.Particle && inter_j.Value == Interaction.Particle)
                        {
                            if (d_j > 0)
                            {
                                force += -0.5f * (pr_i + pr_j) / d_j * SpikyKernelDerivative(p_i, p_j);
                                force += viscosity * (v_j.Value - v_i.Value) / d_j * ViscosityKernel(p_i, p_j);
                            }
                        }

                        if ((inter_i.Value == Interaction.Particle && inter_j.Value == Interaction.Object) ||
                            (inter_i.Value == Interaction.Object && inter_j.Value == Interaction.Particle))
                        {
                            force += -repulsion * ObjectKernel(p_i, p_j);
                            //if (d_j > 0)
                            //{
                            //    force += -3f * (pr_i + pr_j) / d_j * SpikyKernelDerivative(p_i, p_j);
                            //    force += viscosity * (v_j.Value - v_i.Value) / d_j * ViscosityKernel(p_i, p_j);
                            //}
                        }
                    }

                    while (cell_positions.TryGetNextValue(out p_j, ref pos_it) &&
                           hash_map.TryGetNextValue(out j, ref hash_it) &&
                           cell_velocities.TryGetNextValue(out v_j, ref vel_it) && 
                           cell_interactions.TryGetNextValue(out inter_j, ref inter_it))
                    {
                        if (i != j)
                        {
                            float d_j = GetDensity(cell, j);
                            float pr_j = gas_constant * (d_j - ref_density);

                            if (inter_i.Value == Interaction.Particle && inter_j.Value == Interaction.Particle)
                            {
                                if (d_j > 0)
                                {
                                    force += -0.5f * (pr_i + pr_j) / d_j * SpikyKernelDerivative(p_i, p_j);
                                    force += viscosity * (v_j.Value - v_i.Value) / d_j * ViscosityKernel(p_i, p_j);
                                }
                            }

                            if ((inter_i.Value == Interaction.Particle && inter_j.Value == Interaction.Object) ||
                                (inter_i.Value == Interaction.Object && inter_j.Value == Interaction.Particle))
                            {
                                force += -repulsion * ObjectKernel(p_i, p_j);
                                //if (d_j > 0)
                                //{
                                //    force += -3f * (pr_i + pr_j) / d_j * SpikyKernelDerivative(p_i, p_j);
                                //    force += viscosity * (v_j.Value - v_i.Value) / d_j * ViscosityKernel(p_i, p_j);
                                //}
                            }
                        }
                    }
                }

                float max_force = 100f;
                if (Mathf.Abs(force.x) > max_force)
                    force.x = max_force * Mathf.Sign(force.x);
                if (Mathf.Abs(force.y) > max_force)
                    force.y = max_force * Mathf.Sign(force.y);
                if (Mathf.Abs(force.z) > max_force)
                    force.z = max_force * Mathf.Sign(force.z);
                return force;
            }

            public float3 SumForces(int i, Position p_i, Velocity v_i, InteractionType inter_i, float pr_i, int3 cell)
            {
                float3 result = new float3(0f, 0f, 0f);

                for (int dz = -1; dz <= 1; dz++)
                {
                    int cell_index_z = cell.z + dz;

                    for (int dy = -1; dy <= 1; dy++)
                    {
                        int cell_index_y = cell.y + dy;

                        for (int dx = -1; dx <= 1; dx++)
                        {
                            int cell_index_x = cell.x + dx;

                            int current_cell = cell_index_z * hash_constant * hash_constant +
                                cell_index_y * hash_constant +
                                cell_index_x;

                            result += SumForcesCell(i, p_i, v_i, inter_i, pr_i, current_cell);
                        }
                    }
                }

                return result;
            }

            public void Execute(int id)
            {
                if (set_filled[id] == 0)
                    return;

                int cell_hash = set_hash[id];

                int cell_index_x = cell_hash % hash_constant;
                int cell_index_y = cell_hash % (hash_constant * hash_constant);
                cell_index_y /= hash_constant;
                int cell_index_z = cell_hash / (hash_constant * hash_constant);

                int i;
                Position p_i;
                Velocity v_i;
                InteractionType inter_i;

                NativeMultiHashMapIterator<int> pos_it;
                NativeMultiHashMapIterator<int> vel_it;
                NativeMultiHashMapIterator<int> inter_it;
                NativeMultiHashMapIterator<int> hash_it;

                if (cell_positions.TryGetFirstValue(cell_hash, out p_i, out pos_it) &&
                    cell_velocities.TryGetFirstValue(cell_hash, out v_i, out vel_it) &&
                    hash_map.TryGetFirstValue(cell_hash, out i, out hash_it) && 
                    cell_interactions.TryGetFirstValue(cell_hash, out inter_i, out inter_it))
                {
                    int3 cell_index = new int3(
                        cell_index_x,
                        cell_index_y,
                        cell_index_z);

                    float d_i = GetDensity(cell_hash, i);
                    float pr_i = gas_constant * (d_i - ref_density);

                    float3 result = SumForces(i, p_i, v_i, inter_i, pr_i, cell_index);

                    Force force = new Force
                    {
                        Value = result
                    };

                    cell_forces.Add(cell_hash, force);
                    force_hash_map.Add(cell_hash, i);

                    while (cell_positions.TryGetNextValue(out p_i, ref pos_it) &&
                           cell_velocities.TryGetNextValue(out v_i, ref vel_it) &&
                           hash_map.TryGetNextValue(out i, ref hash_it) &&
                           cell_interactions.TryGetNextValue(out inter_i, ref inter_it))
                    {
                        d_i = GetDensity(cell_hash, i);
                        pr_i = gas_constant * (d_i - ref_density);

                        result = SumForces(i, p_i, v_i, inter_i, pr_i,  cell_index);

                        force = new Force
                        {
                            Value = result
                        };

                        cell_forces.Add(cell_hash, force);
                        force_hash_map.Add(cell_hash, i);
                    }
                }

            }
        }

        [Unity.Burst.BurstCompileAttribute]
        struct UnhashForcesJob : IJobParallelFor
        {
            [ReadOnly] public NativeArray<int> hashes;
            [ReadOnly] public NativeMultiHashMap<int, int> hash_map;
            [ReadOnly] public NativeMultiHashMap<int, Force> cell_forces;

            public NativeArray<Force> forces;

            public void Execute(int i)
            {
                NativeMultiHashMapIterator<int> hash_it;
                NativeMultiHashMapIterator<int> force_it;

                int j;
                Force force;
                int cell = hashes[i];

                if (hash_map.TryGetFirstValue(cell, out j, out hash_it) &&
                    cell_forces.TryGetFirstValue(cell, out force, out force_it))
                {

                    if (i == j)
                    {
                        forces[i] = force;
                        return;
                    }

                    while (hash_map.TryGetNextValue(out j, ref hash_it) &&
                           cell_forces.TryGetNextValue(out force, ref force_it))
                    {
                        if (i == j)
                        {
                            forces[i] = force;
                            return;
                        }
                    }
                }

            }
        }

        [Unity.Burst.BurstCompileAttribute]
        struct AddExternalForces : IJobParallelFor
        {
            public NativeArray<Force> forces;

            [ReadOnly] public ComponentDataArray<Position> positions;
            [ReadOnly] public ComponentDataArray<Velocity> velocities;
            [ReadOnly] public float wall_damping;
            [ReadOnly] public float3 system_size;
            [ReadOnly] public float gravity;

            [ReadOnly] public float time;
            [ReadOnly] public float wave_magnitude;
            [ReadOnly] public float wave_period;
            [ReadOnly] public float wave_length;

            public float ComputeRepulsion(float depth, float vel)
            {
                float repulsive_depth = 0.1f;
                float repulsive_force = 10f;

                float rho = (depth / repulsive_depth);

                return repulsive_force * rho - wall_damping * vel;
            }

            public void Execute(int i)
            {
                // Boundary fources
                float3 pos = positions[i].Value;
                float3 force = forces[i].Value;

                if (pos.x < 0f)
                    force.x += ComputeRepulsion(-pos.x, velocities[i].Value.x);
                if (pos.y < 0f)
                    force.y += ComputeRepulsion(-pos.y, velocities[i].Value.y);
                if (pos.z < 0f)
                    force.z += ComputeRepulsion(-pos.z, velocities[i].Value.z);

                if (pos.x > system_size.x)
                    force.x += -ComputeRepulsion(pos.x - system_size.x, velocities[i].Value.x);
                if (pos.y > system_size.y)
                    force.y += -ComputeRepulsion(pos.y - system_size.y, velocities[i].Value.y);
                if (pos.z > system_size.z)
                    force.z += -ComputeRepulsion(pos.z - system_size.z, velocities[i].Value.z);

                force.y += wave_magnitude * Mathf.Sin(2 * Mathf.PI * (pos.x / wave_length - time / wave_period));

                // Gravity
                force.y -= gravity;

                forces[i] = new Force { Value = force };
            }
        }

        [Unity.Burst.BurstCompileAttribute]
        struct GradientDescentUpdateJob : IJobParallelFor
        {
            public ComponentDataArray<Position> positions;

            [ReadOnly] public NativeArray<Force> forces;
            [ReadOnly] public Vector3 size;
            [ReadOnly] public float dt;

            public void Execute(int i)
            {
                float3 pos = positions[i].Value;
                pos += forces[i].Value * dt;

                positions[i] = new Position { Value = pos };
            }
        }

        [Unity.Burst.BurstCompileAttribute]
        struct FinalizeVelocityVerletStepJob : IJobParallelFor
        {
            [ReadOnly] public ComponentDataArray<Force> forces;
            [ReadOnly] public float dt;

            public ComponentDataArray<Velocity> velocities;

            public void Execute(int i)
            {
                velocities[i] = new Velocity
                {
                    Value = velocities[i].Value + 0.5f * dt * forces[i].Value
                };
            }
        }

        [Unity.Burst.BurstCompileAttribute]
        struct CopyForcesJob : IJobParallelFor
        {
            [ReadOnly] public NativeArray<Force> source_forces;

            public ComponentDataArray<Force> dest_forces;

            public void Execute(int i)
            {
                dest_forces[i] = source_forces[i];
            }
        }

        [Unity.Burst.BurstCompileAttribute]
        struct AddBondForcesJob : IJob
        {
            [ReadOnly] public ComponentDataArray<Position> positions;
            [ReadOnly] public ComponentDataArray<Velocity> velocities;
            [ReadOnly] public ComponentDataArray<BondID> ids;
            [ReadOnly] public NativeHashMap<int, BondDistance> bond_distances;

            public ComponentDataArray<Force> forces;

            [ReadOnly] public float bond_strength;

            public float3 ComputeForce(float3 dr, float distance)
            {
                float mag = Mathf.Sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);

                // NOTE(schsam): Optimize.
                return bond_strength / distance / distance * (distance - mag) * dr / mag;
            }

            public void Execute()
            {
                for (int i = 0; i < positions.Length; i++)
                {
                    for (int j = i + 1; j < positions.Length; j++)
                    {
                        float3 dr = positions[i].Value - positions[j].Value;
                        int bid_i = ids[i].Value;
                        int bid_j = ids[j].Value;

                        int hash = bid_i * hash_constant + bid_j;
                        BondDistance bd;

                        if (bond_distances.TryGetValue(hash, out bd))
                        {
                            float3 force = ComputeForce(dr, bd.Value);
                            forces[i] = new Force
                            {
                                Value = forces[i].Value + force
                            };

                            forces[j] = new Force
                            {
                                Value = forces[j].Value - force
                            };
                        }
                    }

                    float damping = 0.9f;
                    forces[i] = new Force
                    {
                        Value = forces[i].Value - damping * velocities[i].Value
                    };

                    /*
                    float max_force = 100f;
                    float3 forcei = forces[i].Value;
                    if (Mathf.Abs(forcei.x) > max_force)
                        forcei.x = max_force * Mathf.Sign(forcei.x);
                    if (Mathf.Abs(forcei.y) > max_force)
                        forcei.y = max_force * Mathf.Sign(forcei.y);
                    if (Mathf.Abs(forcei.z) > max_force)
                        forcei.z = max_force * Mathf.Sign(forcei.z);

                    forces[i] = new Force
                    {
                        Value = forcei
                    };
                    */
                }
            }
        }

        [Unity.Burst.BurstCompileAttribute]
        struct ComputeInteractionCountJob : IJob
        {
            [ReadOnly] public ComponentDataArray<Position> positions;
            [ReadOnly] public NativeMultiHashMap<int, Position> cell_positions;
            [ReadOnly] public NativeMultiHashMap<int, InteractionType> cell_interactions;
            [ReadOnly] public float size_squared;
            [ReadOnly] public float cell_size;

            public NativeArray<int> interaction_count;

            public int CountInteractionsCell(Position p_i, int cell)
            {
                int count = 0;
                Position p_j;
                InteractionType inter_j;

                NativeMultiHashMapIterator<int> pos_it;
                NativeMultiHashMapIterator<int> inter_it;

                if (cell_positions.TryGetFirstValue(cell, out p_j, out pos_it) &&
                    cell_interactions.TryGetFirstValue(cell, out inter_j, out inter_it))
                {
                    if (inter_j.Value == Interaction.Particle)
                    {
                        float3 dr = p_i.Value - p_j.Value;
                        float normSq = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;

                        if (normSq < size_squared)
                            count++;
                    }

                    while (cell_positions.TryGetNextValue(out p_j, ref pos_it) &&
                           cell_interactions.TryGetNextValue(out inter_j, ref inter_it))
                    {
                        if (inter_j.Value == Interaction.Particle)
                        {
                            float3 dr = p_i.Value - p_j.Value;
                            float normSq = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;

                            if (normSq < size_squared)
                                count++;
                        }
                    }
                }

                return count;
            }

            public int CountInteractions(Position p_i, int3 cell)
            {
                int count = 0;
                for (int dz = -1; dz <= 1; dz++)
                {
                    int cell_index_z = cell.z + dz;

                    for (int dy = -1; dy <= 1; dy++)
                    {
                        int cell_index_y = cell.y + dy;

                        for (int dx = -1; dx <= 1; dx++)
                        {
                            int cell_index_x = cell.x + dx;

                            int current_cell = cell_index_z * hash_constant * hash_constant +
                                cell_index_y * hash_constant +
                                cell_index_x;

                            count += CountInteractionsCell(p_i, current_cell);
                        }
                    }
                }

                return count;
            }

            public void Execute()
            {

                for (int i = 0; i < positions.Length; i++)
                {
                    int x_index = (int)Math.Floor(positions[i].Value.x / cell_size) + 10;
                    int y_index = (int)Math.Floor(positions[i].Value.y / cell_size) + 50;
                    int z_index = (int)Math.Floor(positions[i].Value.z / cell_size) + 10;

                    int3 cell = new int3(x_index, y_index, z_index);

                    interaction_count[i] = CountInteractions(positions[i], cell);
                }
            }

        }

        public enum SimulationMode { NVT, NVE, GD };

        ComponentGroup simulated_group;
        ComponentGroup player_group;


        static public Vector3 system_size;
        static public int particle_count = 0;
        static public float cell_size;
        static public SimulationMode mode;
        static public float particle_size;
        static public float gas_constant;
        static public float reference_density;
        static public float gravity;
        static public float viscosity;
        static public GameObject player;
        public static bool clear_simulation = false;

        protected override JobHandle OnUpdate(JobHandle inputDeps)
        {
            float system_volume = system_size.x * system_size.y * system_size.z;
            float step_size = 1f;

            int3 cells_per_side = new int3(
                (int)Math.Ceiling(system_size.x / cell_size),
                (int)Math.Ceiling(system_size.y / cell_size),
                (int)Math.Ceiling(system_size.z / cell_size));
            float cell_volume = cell_size * cell_size * cell_size;

            int num_cells = (int)Math.Ceiling(system_volume / cell_volume);

            var particle_forces =
                new NativeArray<Force>(particle_count, Allocator.TempJob);

            // Map from particle to hashes.
            var particle_hashes =
                new NativeArray<int>(particle_count, Allocator.TempJob);

            // NOTE(schsam): We need the extra hash map because the order of
            // forces and energies in the bucket will not be the same as the
            // order of positions and velocities.

            // Maps from hashes to particle indices.
            var force_hash_map =
                new NativeMultiHashMap<int, int>(particle_count, Allocator.TempJob);
            var hash_map =
                new NativeMultiHashMap<int, int>(particle_count, Allocator.TempJob);
            var density_hash_map =
                new NativeMultiHashMap<int, int>(particle_count, Allocator.TempJob);

            // Set of hashes.
            // We could probably do a better job of picking these sizes.
            int hash_set_size = particle_count * 2;
            var set_hash = new NativeArray<int>(hash_set_size, Allocator.TempJob);
            var set_filled = new NativeArray<int>(hash_set_size, Allocator.TempJob);
            var unique_hash_count = new NativeArray<int>(1, Allocator.TempJob); 

            // All hashed quantities.
            var cell_positions =
                new NativeMultiHashMap<int, Position>(particle_count, Allocator.TempJob);
            var cell_velocities =
                new NativeMultiHashMap<int, Velocity>(particle_count, Allocator.TempJob);
            var cell_forces =
                new NativeMultiHashMap<int, Force>(particle_count, Allocator.TempJob);
            var cell_densities =
                new NativeMultiHashMap<int, float>(particle_count, Allocator.TempJob);
            var cell_interactions =
                new NativeMultiHashMap<int, InteractionType>(particle_count, Allocator.TempJob);

            // Player things
            Player p = player.GetComponent<Player>();

            var player_interactions = new NativeArray<int>(p.positions.Length, Allocator.TempJob, NativeArrayOptions.ClearMemory);

            float dt = step_size * Time.deltaTime;

            JobHandle current_deps = inputDeps;
            JobHandle compute_nose_hoover_zeta_finalize_handle = inputDeps;

            var initialize_hash_set_job = new InitializeHashSet()
            {
                set_hash = set_hash,
                set_filled = set_filled
            };

            var initialize_hash_set_handle =
                initialize_hash_set_job.Schedule(hash_set_size, 64, inputDeps);

            if (mode == SimulationMode.NVE)
            {
                var velocity_verlet_half_step_job = new ComputeVelocityVerletHalfStepJob()
                {
                    forces = simulated_group.GetComponentDataArray<Force>(),
                    dt = dt,
                    size = system_size,
                    velocities = simulated_group.GetComponentDataArray<Velocity>(),
                    positions = simulated_group.GetComponentDataArray<Position>()
                };

                current_deps =
                    velocity_verlet_half_step_job.Schedule(
                        particle_count, 64, current_deps);
            }
            else if (mode == SimulationMode.NVT)
            {
               
            }

            var hash_particle_job = new HashParticlesJob()
            {
                positions = simulated_group.GetComponentDataArray<Position>(),
                velocities = simulated_group.GetComponentDataArray<Velocity>(),
                interactions = simulated_group.GetComponentDataArray<InteractionType>(),
                hash_map = hash_map,
                hash_positions = cell_positions,
                hash_velocities = cell_velocities,
                hash_interactions = cell_interactions,
                hashes = particle_hashes,
                cell_size = cell_size,
                cells_per_side = cells_per_side
            };

            var hash_particle_handle =
                hash_particle_job.Schedule(particle_count, 64, current_deps);

            var hash_deps_handle = JobHandle.CombineDependencies(hash_particle_handle, initialize_hash_set_handle);

            var build_hash_set_job = new BuildHashSetJob()
            {
                hashes = particle_hashes,
                set_hash = set_hash,
                set_filled = set_filled,
                unique_hash_count = unique_hash_count
            };

            //Debug.Log(unique_hash_count);
            var build_hash_set_handle = build_hash_set_job.Schedule(hash_deps_handle);

            var compute_density_job = new ComputeDensityJob()
            {
                hash_map = hash_map,
                cell_positions = cell_positions,
                size = particle_size,
                size2 = particle_size * particle_size,
                isize9 = Mathf.Pow(particle_size, -9f),
                isize6 = Mathf.Pow(particle_size, -6f),
                cells_per_side = cells_per_side,
                num_cells = num_cells,
                system_size = system_size,
                set_hash = set_hash,
                set_filled = set_filled,
                density_hash_map = density_hash_map,
                cell_density = cell_densities,
                cell_interactions = cell_interactions
            };

            var compute_density_handle =
                compute_density_job.Schedule(hash_set_size, 64, build_hash_set_handle);

            var compute_forces_job = new ComputeForcesJob()
            {
                hash_map = hash_map,
                cell_positions = cell_positions,
                cell_velocities = cell_velocities,
                cell_interactions = cell_interactions,
                size = particle_size,
                size2 = particle_size * particle_size,
                isize9 = Mathf.Pow(particle_size, -9f),
                isize6 = Mathf.Pow(particle_size, -6f),
                cells_per_side = cells_per_side,
                num_cells = num_cells,
                cell_forces = cell_forces,
                force_hash_map = force_hash_map,
                system_size = system_size,
                set_hash = set_hash,
                set_filled = set_filled,
                density_hash_map = density_hash_map,
                cell_densities = cell_densities,
                gas_constant = gas_constant,
                ref_density = reference_density,
                viscosity = viscosity,
                repulsion = 20f
            };

            var compute_forces_handle =
                compute_forces_job.Schedule(hash_set_size, 64, compute_density_handle);

            var unhash_forces_job = new UnhashForcesJob()
            {
                hashes = particle_hashes,
                hash_map = force_hash_map,
                cell_forces = cell_forces,
                forces = particle_forces,
            };

            var unhash_forces_handle =
                unhash_forces_job.Schedule(particle_count, 64, compute_forces_handle);

            var add_external_forces_job = new AddExternalForces()
            {
                positions = simulated_group.GetComponentDataArray<Position>(),
                velocities = simulated_group.GetComponentDataArray<Velocity>(),
                forces = particle_forces,
                system_size = system_size,
                gravity = gravity,
                wall_damping = 0.7f,
                wave_period = 8f,
                wave_magnitude = p.wave_magnitude,
                wave_length = 15f,
                time = Time.time
            };

            var add_external_forces_handle =
                add_external_forces_job.Schedule(particle_count, 64, unhash_forces_handle);

            var copy_forces_job = new CopyForcesJob()
            {
                source_forces = particle_forces,
                dest_forces = simulated_group.GetComponentDataArray<Force>()
            };

            var copy_forces_handle = copy_forces_job.Schedule(particle_count, 64, add_external_forces_handle);

            var add_bond_forces_job = new AddBondForcesJob()
            {
                positions = player_group.GetComponentDataArray<Position>(),
                velocities = player_group.GetComponentDataArray<Velocity>(),
                ids = player_group.GetComponentDataArray<BondID>(),
                bond_distances = p.bonds,
                forces = player_group.GetComponentDataArray<Force>(),
                bond_strength = 50f
            };
            
            var add_bond_forces_handle = add_bond_forces_job.Schedule(copy_forces_handle);

            var count_interactions_job = new ComputeInteractionCountJob()
            {
                positions = player_group.GetComponentDataArray<Position>(),
                cell_positions = cell_positions,
                cell_interactions = cell_interactions,
                cell_size = cell_size,
                size_squared = particle_size * particle_size,
                interaction_count = player_interactions
            };

            JobHandle count_interactions_handle = count_interactions_job.Schedule(add_bond_forces_handle);

            JobHandle step_finalize_handle = add_bond_forces_handle;

            if (mode == SimulationMode.NVE)
            {
                var velocity_verlet_finalize_job = new FinalizeVelocityVerletStepJob()
                {
                    dt = dt,
                    velocities = simulated_group.GetComponentDataArray<Velocity>(),
                    forces = simulated_group.GetComponentDataArray<Force>()
                };

                step_finalize_handle =
                    velocity_verlet_finalize_job.Schedule(particle_count, 64, step_finalize_handle);

            }
            else if (mode == SimulationMode.NVT)
            {

            }
            else if (mode == SimulationMode.GD)
            {
                var position_job = new GradientDescentUpdateJob()
                {
                    dt = dt,
                    positions = simulated_group.GetComponentDataArray<Position>(),
                    forces = particle_forces,
                    size = system_size
                };

                step_finalize_handle =
                    position_job.Schedule(particle_count, 64, step_finalize_handle);
            }

            step_finalize_handle.Complete();
            count_interactions_handle.Complete();


            if (!p.has_lost)
            {
                int num_no_interactions = 0;
                for (int i = 0; i < player_interactions.Length; i++)
                {
                    if (player_interactions[i] == 0)
                        num_no_interactions++;
                }

                if (num_no_interactions < 3)
                {
                    if (!p.submerged)
                        p.submerged_time = Time.time;

                    float sdt = Time.time - p.submerged_time;

                    if (sdt > 1f)
                    {
                        p.submerged_object.GetComponent<MeshRenderer>().enabled = true;

                        Transform t = p.submerged_object.GetComponent<Transform>();
                        t.localScale = new Vector3(10f * (1f - sdt / p.max_submerged_time), t.localScale.y, t.localScale.z);
                        t.localPosition = new Vector3(t.localScale.x / 2f - 0.5f, t.localPosition.y, t.localPosition.z);

                        var c = p.submerged_text_object.GetComponent<TMPro.TMP_Text>().color;
                        c.a = 1f;
                        p.submerged_text_object.GetComponent<TMPro.TMP_Text>().color = c;
                    }

                    if (sdt > p.max_submerged_time)
                    {
                        p.has_lost = true;

                        var c = p.lost_text_object.GetComponent<TMPro.TMP_Text>().color;
                        c.a = 1f;
                        p.lost_text_object.GetComponent<TMPro.TMP_Text>().color = c;

                        c = p.goal_1_object.GetComponent<TMPro.TMP_Text>().color;
                        c.a = 0f;
                        p.goal_1_object.GetComponent<TMPro.TMP_Text>().color = c;

                        c = p.goal_2_object.GetComponent<TMPro.TMP_Text>().color;
                        c.a = 0f;
                        p.goal_2_object.GetComponent<TMPro.TMP_Text>().color = c;

                        p.submerged_object.GetComponent<MeshRenderer>().enabled = false;

                        c = p.restart_text_object.GetComponent<TMPro.TMP_Text>().color;
                        c.a = 1f;
                        p.restart_text_object.GetComponent<TMPro.TMP_Text>().color = c;
                    }

                    p.submerged = true;
                }
                else
                {
                    p.submerged = false;

                    p.submerged_object.GetComponent<MeshRenderer>().enabled = false;

                    var c = p.submerged_text_object.GetComponent<TMPro.TMP_Text>().color;
                    c.a = 0f;
                    p.submerged_text_object.GetComponent<TMPro.TMP_Text>().color = c;
                }
            }

            unique_hash_count.Dispose();
            set_hash.Dispose();
            set_filled.Dispose();
            cell_velocities.Dispose();
            force_hash_map.Dispose();
            particle_forces.Dispose();
            particle_hashes.Dispose();
            hash_map.Dispose();
            cell_positions.Dispose();
            cell_forces.Dispose();
            cell_densities.Dispose();
            density_hash_map.Dispose();
            cell_interactions.Dispose();
            player_interactions.Dispose();

            return hash_particle_handle;
        }

        protected override void OnCreateManager(int capacity)
        {
            simulated_group = GetComponentGroup(
                typeof(Position),
                typeof(Velocity),
                typeof(Force),
                typeof(InteractionType));

            player_group = GetComponentGroup(
                typeof(Position),
                typeof(Velocity),
                typeof(Force),
                typeof(BondID));
        }

    }
}