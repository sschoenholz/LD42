using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraControl : MonoBehaviour {

    public Transform target;
    public Transform simulation;

    // Use this for initialization
    void Start()
    {
        target.parent = simulation;
        transform.parent = target;
    }

    // Update is called once per frame
    void Update()
    {
        transform.LookAt(target);
    }
}
