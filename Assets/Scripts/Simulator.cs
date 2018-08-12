using System;
using Unity.Entities;
using Unity.Rendering;
using UnityEngine;
using Unity.Transforms;
using Unity.Mathematics;

namespace ParticleSimulator
{
    public struct Particle : IComponentData { }

    public enum Interaction { Particle, Object };
    public struct InteractionType : IComponentData { public Interaction Value; }

    public struct Velocity : IComponentData { public float3 Value; }

    public struct Force : IComponentData { public float3 Value; }

    public struct Color : IComponentData { public float4 Value; }

    public struct BondID : IComponentData { public int Value; }

    public struct BondDistance : IComponentData { public float Value; }

    public class Simulator : MonoBehaviour
    {
        [SerializeField]
        public int particle_count = 10000;

        [SerializeField]
        public float particle_radius = 1.0f;

        [SerializeField]
        public float cell_size = 2.0f;

        [SerializeField]
        public float gas_constant = 1f;
        
        [SerializeField]
        public float reference_density = 20f;

        [SerializeField]
        public float viscosity = 0.05f;

        [SerializeField]
        public float gravity = 0.1f;

        [SerializeField]
        public Vector3 size = new Vector3(10f, 10f, 10f);

        [SerializeField]
        public GameObject particle_object;

        [SerializeField]
        public GameObject particle_object_s2;

        [SerializeField]
        public SimulateFrame.SimulationMode mode;

        [SerializeField]
        public float temperature;

        [SerializeField]
        public float time_between_pour;

        [SerializeField]
        public GameObject player; 

        EntityArchetype particle_archetype;

        MeshInstanceRenderer particle_renderer;

        bool pouring;

        float last_pour_time;

        // Use this for initialization
        void Start()
        {
            var entity_manager = World.Active.GetOrCreateManager<EntityManager>();

            particle_archetype = entity_manager.CreateArchetype(
                typeof(Particle), typeof(Position), typeof(Velocity), typeof(Force), typeof(TransformMatrix), typeof(Color), typeof(InteractionType));

            var particle_look = GameObject.Instantiate(particle_object);
            particle_renderer = particle_look.GetComponent<MeshInstanceRendererComponent>().Value;

            particle_renderer.receiveShadows = true;
            particle_renderer.castShadows = UnityEngine.Rendering.ShadowCastingMode.On;

            float velocity_rescale = Mathf.Sqrt(3f * temperature);
            for (int i = 0; i < particle_count; i++)
            {
                Entity particle = entity_manager.CreateEntity(particle_archetype);

                float3 position = new float3(UnityEngine.Random.Range(0f, size.x),
                                       UnityEngine.Random.Range(0f, 2.0f),
                                       UnityEngine.Random.Range(0f, size.z));

                entity_manager.SetComponentData(particle, new Position
                {
                    Value = new float3(position)
                });

                entity_manager.SetComponentData(particle, new Velocity
                {
                    Value = velocity_rescale * (new float3(UnityEngine.Random.Range(-1f, 1f),
                                       UnityEngine.Random.Range(-0.1f, 0.1f),
                                       UnityEngine.Random.Range(-1f, 1f)))
                });

                entity_manager.SetComponentData(particle, new Force
                {
                    Value = new float3(0f, 0f, 0f)
                });

                entity_manager.SetComponentData(particle, new Color
                {
                    Value = new float4(0f, 1f, 1f, 1f)
                });

                entity_manager.SetComponentData(particle, new InteractionType
                {
                    Value = Interaction.Particle
                });

                entity_manager.AddSharedComponentData(particle, particle_renderer);
            }

            GameObject.Destroy(particle_look);

            SimulateFrame.system_size = size;
            SimulateFrame.particle_count += particle_count;
            SimulateFrame.cell_size = cell_size;
            SimulateFrame.mode = mode;
            SimulateFrame.particle_size = particle_radius;
            SimulateFrame.reference_density = reference_density;
            SimulateFrame.gravity = gravity;
            SimulateFrame.gas_constant = gas_constant;
            SimulateFrame.viscosity = viscosity;
            SimulateFrame.player = player;

            Application.targetFrameRate = 60;

            current_mode = mode;

            vertical_rotation = 0f;
            last_pour_time = 0;

            water_position = new Vector3(5f, 5f, 5f);

            pouring = false;
        }

        SimulateFrame.SimulationMode current_mode;
        float vertical_rotation;

        Vector3 water_position;

        [SerializeField]
        public float water_speed;

        private void PourWater(Vector3 pos, bool type)
        {
            var entity_manager = World.Active.GetOrCreateManager<EntityManager>();

            particle_archetype = entity_manager.CreateArchetype(
                typeof(Particle), typeof(Position), typeof(Velocity), typeof(Force), typeof(TransformMatrix), typeof(Color), typeof(InteractionType));

            GameObject particle_look;
            if (type)
                particle_look = GameObject.Instantiate(particle_object_s2);
            else
                particle_look = GameObject.Instantiate(particle_object);

            particle_renderer = particle_look.GetComponent<MeshInstanceRendererComponent>().Value;

            particle_renderer.receiveShadows = false;
            particle_renderer.castShadows = UnityEngine.Rendering.ShadowCastingMode.Off;

            float velocity_rescale = Mathf.Sqrt(3f * temperature);

            Player p = player.GetComponent<Player>();

            for (int i = 0; i < p.water_per_pour; i++)
            {
                Entity particle = entity_manager.CreateEntity(particle_archetype);

                float3 position = new float3(UnityEngine.Random.Range(pos.x - 0.5f, pos.x + 0.5f),
                                        UnityEngine.Random.Range(20f, 22f),
                                       UnityEngine.Random.Range(pos.z - 0.5f, pos.z + 0.5f));
                entity_manager.SetComponentData(particle, new Position
                {
                    Value = new float3(position)
                });


                entity_manager.SetComponentData(particle, new Velocity
                {
                    Value = velocity_rescale * (new float3(0f, -20f, 0f))
                });

                entity_manager.SetComponentData(particle, new Force
                {
                    Value = new float3(0f, 0f, 0f)

                });

                entity_manager.SetComponentData(particle, new Color
                {
                    Value = new float4(0f, 1f, 1f, 1f)
                });

                entity_manager.SetComponentData(particle, new InteractionType
                {
                    Value = Interaction.Particle
                });

                entity_manager.AddSharedComponentData(particle, particle_renderer);
            }
            SimulateFrame.particle_count += p.water_per_pour;
            GameObject.Destroy(particle_look);

        }

        // Update is called once per frame
        void Update()
        {
            Player p = player.GetComponent<Player>();

            if (p.has_won || p.has_lost)
                return;

            if (Input.GetButtonDown("Pour"))
            {
                pouring = true;
            }


           if ((pouring || p.water_leaking > 0) && Time.time - last_pour_time > time_between_pour)
            {
                last_pour_time = Time.time;
                PourWater(new Vector3(10f, 0f, 10f), false);

                if (p.water_leaking > 1)
                {
                    Vector3 dr = player.GetComponent<Transform>().position - water_position;
                    dr.Normalize();

                    water_position += water_speed * dr * Time.deltaTime;

                    PourWater(water_position, true);
                }
            }

        }

        private void OnDestroy()
        {
        }
    }

}
