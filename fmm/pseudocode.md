Function FMM_Gravitational_Forces(particles, theta)
    Input: 
        particles: list of N particles with (position, mass)
        theta: opening angle parameter (controls approximation accuracy)
    Output:
        forces: list of 3D force vectors on each particle

    // Step 1: Build Octree
    root = BuildOctree(particles)

    // Step 2: Upward Pass - Multipole Expansion
    Function ComputeMultipoles(node)
        if node is leaf:
            ComputeMultipoleForLeaf(node)
        else:
            for child in node.children:
                ComputeMultipoles(child)
            AggregateChildMultipolesToParent(node)
    end
    ComputeMultipoles(root)

    // Step 3: Downward Pass - Local Expansions
    Function ComputeLocalExpansions(node, parent_interactions)
        interaction_list = GetInteractionList(node, theta)
        for src_node in interaction_list:
            TranslateMultipoleToLocal(node, src_node)
        
        if node is not leaf:
            for child in node.children:
                ComputeLocalExpansions(child, node.local_expansion)
    end
    ComputeLocalExpansions(root, null)

    // Step 4: Evaluate Forces
    forces = []
    for each leaf_node in OctreeLeaves(root):
        for particle in leaf_node.particles:
            // Near-field: direct summation
            force_direct = 0
            for neighbor in GetNeighboringLeafNodes(leaf_node):
                for p2 in neighbor.particles:
                    if p2 ≠ particle:
                        force_direct += ComputeDirectForce(particle, p2)

            // Far-field: use local expansion
            force_far = EvaluateLocalExpansion(leaf_node, particle.position)

            total_force = force_direct + force_far
            forces.append(total_force)

    return forces
End Function
