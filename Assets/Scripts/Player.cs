using System;
using Unity.Entities;
using UnityEngine;
using Unity.Transforms;
using Unity.Mathematics;
using Unity.Collections;
using System.Collections;
using UnityEngine.SceneManagement;

namespace ParticleSimulator
{


    public class Player : MonoBehaviour {
        
        public const int hash_constant = 1997;

        public GameObject[] debug_objects;
        public Entity[] physical;
        public float3[] positions;
        public NativeHashMap<int, BondDistance> bonds;

        public bool submerged;

        [SerializeField]
        public GameObject debug_prototype;

        public float3 center;

        public EntityArchetype player_archetype;

        public bool thrusters;

        [SerializeField]
        public float thrust_size;
        public float fuel;

        [SerializeField]
        public float max_fuel;
        [SerializeField]
        public float fuel_burn_rate;
        [SerializeField]
        public float fuel_recharge_rate;

        public float start_time;

        [SerializeField]
        public float message_display_time;
        [SerializeField]
        public float message_fade_out_lag;
        [SerializeField]
        public float message_fade_in_lag;

        [SerializeField]
        public float press_key_message_time;
        [SerializeField]
        public GameObject press_key_object;

        [SerializeField]
        public GameObject title_object;

        [SerializeField]
        public GameObject story_1_object;
        [SerializeField]
        public GameObject story_2_object;
        [SerializeField]
        public GameObject story_3_object;
        [SerializeField]
        public GameObject story_4_object;
        [SerializeField]
        public GameObject story_5_object;
        [SerializeField]
        public GameObject story_6_object;
        [SerializeField]
        public GameObject story_7_object;
        [SerializeField]
        public GameObject story_8_object;

        [SerializeField]
        public GameObject goal_1_object;
        [SerializeField]
        public GameObject goal_2_object;

        [SerializeField]
        public GameObject victory_1_object;
        [SerializeField]
        public GameObject victory_2_object;

        [SerializeField]
        public GameObject uh_object;
        [SerializeField]
        public GameObject oh_object;

        [SerializeField]
        public GameObject submerged_text_object;
        [SerializeField]
        public GameObject submerged_object;

        [SerializeField]
        public GameObject lost_text_object;
        [SerializeField]
        public GameObject restart_text_object;


        public float submerged_time;
        [SerializeField]
        public float max_submerged_time;

        public float wave_magnitude;

        int current_story;

        bool press_key_displayed;
        bool key_pressed;

        bool engines_work;

        public int water_leaking;
        public int water_per_pour;

        bool finished_tutorial;
        public bool has_lost;
        public bool has_won;

        // Use this for initialization
        void Start()
        {
            Transform t = GetComponent<Transform>();

            Vector3 scale = t.localScale;

            var entity_manager = World.Active.GetOrCreateManager<EntityManager>();

            player_archetype = entity_manager.CreateArchetype(
                typeof(Position), typeof(Velocity), typeof(Force), typeof(BondID), typeof(InteractionType));

            float xscale = 1.1f * scale.x / 2f;
            float yscale = scale.y / 2f;
            float zscale = 1.1f * scale.z / 2f;

            positions = new float3[] {
                new float3(-xscale , -yscale , -zscale ),

                new float3(-xscale , -yscale ,  0f),
                new float3(-xscale , -yscale ,  zscale ),

                new float3(0f, -yscale ,  -zscale ),
                new float3( xscale , -yscale , -zscale ),

                new float3(0f, -yscale , zscale ),
                new float3(xscale , -yscale , 0f),

                new float3( xscale , -yscale ,  zscale ),
                new float3( 0f , -yscale ,  0f ),

                new float3(-xscale , yscale , -zscale ),

                new float3(-xscale , yscale ,  0f),
                new float3(-xscale , yscale ,  zscale ),

                new float3(0f, yscale ,  -zscale ),
                new float3( xscale , yscale , -zscale ),

                new float3(0f, yscale , zscale ),
                new float3(xscale , yscale , 0f),

                new float3( xscale , yscale ,  zscale ),
                new float3( 0f , yscale ,  0f ),
            };

            if (debug_prototype)
            {
                debug_objects = new GameObject[positions.Length];
                for (int i = 0; i < positions.Length; i++)
                {
                    debug_objects[i] = Instantiate(debug_prototype);
                }
            }

            physical = new Entity[positions.Length];

            for (int i = 0; i < positions.Length; i++)
            {
                physical[i] = entity_manager.CreateEntity(player_archetype);

                entity_manager.SetComponentData(physical[i], new Position
                {
                    Value = positions[i] + center
                });

                entity_manager.SetComponentData(physical[i], new Velocity
                {
                    Value = new float3(0f, 0f, 0f)
                });

                entity_manager.SetComponentData(physical[i], new Force
                {
                    Value = new float3(0f, 0f, 0f)
                });

                entity_manager.SetComponentData(physical[i], new BondID
                {
                    Value = i
                });

                entity_manager.SetComponentData(physical[i], new InteractionType
                {
                    Value = Interaction.Object
                });
            }

            SimulateFrame.particle_count += positions.Length;

            bonds = new NativeHashMap<int, BondDistance>(
                positions.Length * positions.Length, Allocator.Persistent);
            for (int i = 0; i < positions.Length; i++)
            {
                for (int j = 0; j < positions.Length; j++)
                {
                    if (i != j)
                    {
                        int hash = i * hash_constant + j;
                        Vector3 dr = positions[i] - positions[j];
                        if (!bonds.TryAdd(hash, new BondDistance
                        {
                            Value = Vector3.Magnitude(dr)
                        }))
                            Debug.Log("Warning! Could not add bond to bonds list.");
                    }
                }
            }

            engines_work = true;
            water_leaking = 0;
            wave_magnitude = 1.4f;

            finished_tutorial = false;

            press_key_displayed = false;
            key_pressed = false;

            has_lost = false;

            UnityEngine.Color c = press_key_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            press_key_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = story_1_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            story_1_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = story_2_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            story_2_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = story_3_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            story_3_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = story_4_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            story_4_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = story_4_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            story_4_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = story_5_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            story_5_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = story_6_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            story_6_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = story_7_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            story_7_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = story_8_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            story_8_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = goal_1_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            goal_1_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = goal_2_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            goal_2_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = victory_1_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            victory_1_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = victory_2_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            victory_2_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = uh_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            uh_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = oh_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            oh_object.GetComponent<TMPro.TMP_Text>().color = c;

            submerged_object.GetComponent<MeshRenderer>().enabled = false;
            
            c = submerged_text_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            submerged_text_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = lost_text_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            lost_text_object.GetComponent<TMPro.TMP_Text>().color = c;

            c = restart_text_object.GetComponent<TMPro.TMP_Text>().color;
            c.a = 0f;
            restart_text_object.GetComponent<TMPro.TMP_Text>().color = c;

            start_time = Time.time;

            water_per_pour = 2;
            tmp_blah = false;
        }

        private void OnDestroy()
        {
            bonds.Dispose();
        }

        IEnumerator FadeIn(GameObject obj)
        {
            Debug.Log("Starting Fade In");

            float start_fade = Time.time;
            var color = obj.GetComponent<TMPro.TMP_Text>().color;

            float dt = Time.time - start_fade;
            while (dt < message_fade_in_lag)
            {
                dt = Time.time - start_fade;
                yield return null;
            }

            Debug.Log("Passed Initial Fade In Wait");

            start_fade = Time.time;
            dt = Time.time - start_fade;
            while (dt < message_display_time)
            {
                color.a = dt / message_display_time;
                obj.GetComponent<TMPro.TMP_Text>().color = color;
                yield return null;
                dt = Time.time - start_fade;
            }

            Debug.Log("Finishing Fade In");
            color.a = 1f;
            obj.GetComponent<TMPro.TMP_Text>().color = color;
        }

        IEnumerator FadeOut(GameObject obj)
        {
            Debug.Log("Starting Fade Out");

            float start_fade = Time.time;

            var color = obj.GetComponent<TMPro.TMP_Text>().color;
            float initial_a = color.a;

            float dt = Time.time - start_fade;
            while (dt < message_fade_out_lag)
            {
                dt = Time.time - start_fade;
                yield return null;
            }

            Debug.Log("Passed Initial Fade Out Wait");

            start_fade = Time.time;
            dt = Time.time - start_fade;
            while (dt < message_display_time)
            {
                color.a = initial_a * (1f - dt / message_display_time);
                obj.GetComponent<TMPro.TMP_Text>().color = color;
                yield return null;
                dt = Time.time - start_fade;
            }

            Debug.Log("Finishing Fade Out");

            color.a = 0f;
            obj.GetComponent<TMPro.TMP_Text>().color = color;
        }

        bool tmp_blah;

        void HandleStory()
        {
            float ctime = Time.time;

            if ((!key_pressed) && 
                (!press_key_displayed) && 
                ctime - start_time > press_key_message_time)
            {
                StartCoroutine("FadeIn", press_key_object);
                press_key_displayed = true;
            }

            if ((!key_pressed) && (Mathf.Abs(Input.GetAxis("Horizontal")) > 0.1f ||
                Mathf.Abs(Input.GetAxis("Vertical")) > 0.1f))
            {
                StopCoroutine("FadeIn");
                StartCoroutine("FadeOut", press_key_object);
                StartCoroutine("FadeOut", title_object);
                StartCoroutine("FadeIn", story_1_object);

                key_pressed = true;
                start_time = ctime;
                current_story = 1;
            }


            if (water_leaking == 0 && key_pressed)
            {
                if (ctime - start_time > message_fade_in_lag + message_display_time)
                {
                    if (current_story == 1)
                    {
                        StartCoroutine("FadeOut", story_1_object);
                        StartCoroutine("FadeIn", story_2_object);
                        start_time = ctime;
                    }

                    if (current_story == 3)
                    {
                        StartCoroutine("FadeOut", story_2_object);
                        StartCoroutine("FadeOut", story_3_object);

                        StartCoroutine("FadeIn", story_4_object);
                        start_time = ctime;
                    }

                    if (current_story == 5)
                    {
                        StartCoroutine("FadeOut", story_4_object);
                        StartCoroutine("FadeOut", story_5_object);

                        engines_work = false;
                        start_time = ctime;
                    }

                    current_story++;
                }

                if (current_story == 2 && ctime - start_time > press_key_message_time)
                {
                    StartCoroutine("FadeIn", story_3_object);
                    start_time = ctime;
                    current_story++;
                }

                if (current_story == 4 && ctime - start_time > press_key_message_time)
                {
                    StartCoroutine("FadeIn", story_5_object);
                    start_time = ctime;
                    current_story++;
                }
            }

            if (!engines_work)
            {
                if (current_story == 8 && ctime - start_time > message_fade_in_lag + message_display_time + 1)
                {
                    StartCoroutine("FadeOut", uh_object);
                    StartCoroutine("FadeOut", oh_object);
                    StartCoroutine("FadeIn", story_6_object);

                    start_time = ctime;
                    current_story++;
                }

                if (current_story == 7 && ctime - start_time > message_fade_in_lag - 1)
                {
                    water_leaking = 1;
                    wave_magnitude = 0.25f;
                    current_story++;
                }

                if (current_story == 7 && ctime - start_time > 1)
                {
                    StartCoroutine("FadeIn", oh_object);
                }

                if (current_story == 6 && ctime - start_time > press_key_message_time)
                {
                    StartCoroutine("FadeIn", uh_object);
                    start_time = ctime;
                    current_story++;
                }
            }

            if (water_leaking == 1)
            {
                if (current_story == 9 && ctime - start_time > 45f - message_fade_in_lag + 1f) {
                    StartCoroutine("FadeIn", story_7_object);
                    StartCoroutine("FadeIn", story_8_object);
                    current_story++;
                }

                if (ctime - start_time > 45f)
                {
                    water_leaking = 2;
                    wave_magnitude = 0.2f;
                    start_time = ctime;
                    current_story++;
                }
            }

            if (water_leaking == 2)
            {
                if (current_story == 11 && ctime - start_time > 15f)
                {
                    StartCoroutine("FadeOut", story_7_object);
                    StartCoroutine("FadeOut", story_8_object);
                    current_story++;
                }

                if (ctime - start_time > 45f)
                {
                    StartCoroutine("FadeIn", goal_1_object);

                    water_leaking = 3;
                    water_per_pour = 4;
                }

            }

            if (water_leaking == 3)
            {
                if ((!tmp_blah) && ctime - start_time > 65f)
                {
                    StartCoroutine("FadeIn", goal_2_object);
                    tmp_blah = true;
                }

            }
        }

        // Update is called once per frame
        void Update() {
            if (Input.GetKeyDown(KeyCode.Escape))
            {
                if (!key_pressed)
                {
                    Application.Quit();
                }
                else
                {
                    var manager = World.Active.GetOrCreateManager<EntityManager>();
                    var entities = manager.GetAllEntities();

                    foreach (var e in entities)
                        manager.DestroyEntity(e);

                    entities.Dispose();
                    SceneManager.LoadScene(SceneManager.GetActiveScene().buildIndex);
                    SimulateFrame.particle_count = 0;
                }
            }

            if (has_won)
                return;

            if (has_lost)
                return;

            HandleStory();

            // Handle input
            float3 force = new float3(
                13.0f * Input.GetAxis("Horizontal"),
                0f,
                13.0f * Input.GetAxis("Vertical"));

            float gravity = -3.15f;
            float scale = 0.85f;
            float lift = Mathf.Sqrt(force.x * force.x + force.z * force.z) * scale / 15f / Mathf.Sqrt(2f);
            force.y += gravity + scale * lift;

            // Adjust the physical particles
            var entity_manager = World.Active.GetOrCreateManager<EntityManager>();

            for (int i = 0; i < positions.Length; i++)
            {
                positions[i] = entity_manager.GetComponentData<Position>(physical[i]).Value;
                float3 f = entity_manager.GetComponentData<Force>(physical[i]).Value;

                entity_manager.SetComponentData<Force>(physical[i], new Force
                {
                    Value = f + force
                });

                if (positions[i].y >= 14.2f)
                {
                    has_won = true;

                    var c = goal_1_object.GetComponent<TMPro.TMP_Text>().color;
                    c.a = 0f;
                    goal_1_object.GetComponent<TMPro.TMP_Text>().color = c;

                    c = goal_2_object.GetComponent<TMPro.TMP_Text>().color;
                    c.a = 0f;
                    goal_2_object.GetComponent<TMPro.TMP_Text>().color = c;

                    c = victory_1_object.GetComponent<TMPro.TMP_Text>().color;
                    c.a = 1f;
                    victory_1_object.GetComponent<TMPro.TMP_Text>().color = c;

                    c = victory_2_object.GetComponent<TMPro.TMP_Text>().color;
                    c.a = 1f;
                    victory_2_object.GetComponent<TMPro.TMP_Text>().color = c;

                    c = restart_text_object.GetComponent<TMPro.TMP_Text>().color;
                    c.a = 1f;
                    restart_text_object.GetComponent<TMPro.TMP_Text>().color = c;
                }

            }
            Transform t = GetComponent<Transform>();

            float3 average_position = new float3(0f, 0f, 0f);

            for (int i = 0; i < positions.Length; i++)
            {
                if (debug_prototype)
                   debug_objects[i].GetComponent<Transform>().position = positions[i];
                average_position += positions[i] / positions.Length;
            }

            t.position = average_position;

            Vector3 dr1 = Vector3.Normalize(positions[0] - positions[2]);
            Vector3 dr2 = Vector3.Normalize(positions[0] - positions[4]);
            
            Vector3 normal = Vector3.Normalize(Vector3.Cross(dr1, dr2));

            t.rotation = Quaternion.LookRotation(dr1, normal);

            Debug.Log(Time.time - start_time);
            Debug.Log(average_position);

        }
    }
}
