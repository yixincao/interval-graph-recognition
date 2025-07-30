package graph.algorithms;

import java.util.*;

import graph.Graph;
import util.*;
import static util.Permutation.*;

/**
 * 
 * @author Yixin Cao and Chenxi Liu (May 2025)
 *
 *         Implementation of major graph search algorithms.
 * 
 *         By default, an ordering σ is from [0..n-1] → V(G); hence, σ[i] = v
 *         meaning that vertex v is the ith visited (counted from 0).
 * 
 *         In particular, the arguments and return values use this ordering.
 *         The only exception is {@code Permutation.inversePermutation}, which
 *         calculates the inverse (not to be confused with the reverse; see
 *         comments in {@code Permutation}) of an ordering.
 *
 *         This is differenet from the default meaning of orderings in the paper
 *         [C21], where all orderings are from V(G) → [1..n] (not [0..n-1]).
 *
 *         References: LBFS [RTL76]; LBFS+ [RTL76]; LBFS↑ [LW14, C21]; MCS
 *         [TY84].
 *
 *         [RTL76] Rose, Donald J., Tarjan, Robert Endre and Lueker, George S..
 *         Algorithmic aspects of vertex elimination on graphs. 
 *         SIAM J. Comput., 5(2):266--283, 1976. http://doi.org/10.1137/0205021
 *         [TY84] Tarjan, Robert Endre and Yannakakis, Mihalis. 
 *         Simple linear-time algorithms to test chordality of graphs, test acyclicity of hypergraphs, and selectively reduce acyclic hypergraphs. 
 *         SIAM J. Comput., 13(3):566--579, 1984. http://doi.org/10.1137/0213035
 *         [TY85] Addendum in the same journal, 14(1):254-255, 1985. http://doi.org/10.1137/0214020
 *         [LW14] Li, Peng and Wu, Yaokun. A four-sweep LBFS recognition algorithm for interval graphs. 
 *         Discrete Math. Theor. Comput. Sci., 16(3):23--50, 2014. http://doi.org/10.46298/dmtcs.2094
 *         [C21] Cao, Yixin, Recognizing (Unit) Interval Graphs by Zigzag Graph Searches, 
 *         Proceedings of the 4th SIAM Symposium on Simplicity in Algorithms (SOSA), 2021, pages 92--106.
 *         http://doi.org/10.1137/1.9781611976496.11
 *         [S91] Simon, Klaus. A new simple linear algorithm to recognize interval graphs. 
 *         Computational Geometry - Methods, Algorithms and Applications, International Workshop on Computational Geometry CG'91, pages 289--308, 1991.
 * 
 */
public class GraphSearch {
    public static boolean DEBUG = true; // to be removed.

    // exposed for JUnit test.
    public static int[][] sortAdjacencyListsTest(int[][] adj, int[] permutation) {
        return sortAdjacencyLists(adj, permutation);
    }
    public static int[] upOrderingTest(int[][] originalLists, int[] previousOrder) {
        return upOrdering(originalLists, previousOrder);
    }

    /**
     *
     * @param g: the input graph
     * @return 
     */
    public static int[] bfs(Graph g) {
        // not implemented.
        return null;
    }

    /**
     * BFS starting from a given {@code source}, avoiding vertices marked as forbidden.
     *
     * {@code parent[i]} is always -1 when {@code forbidden[i] = true}.
     * 
     * @param g: the input graph
     * @param forbidden: whether a vertex is forbidden (when true)
     * @param source: the starting vertex (which may or may not be forbidden)
     *
     * @return the parent of each vertex in the partial BFS tree rooted at {@code source}
     */
    public static int[] bfsWithForbiddenVertices(Graph g, boolean[] forbidden, int source) {
        return bfsWithForbiddenVertices(g.adj, forbidden, source);
    }
    public static int[] bfsWithForbiddenVertices(int[][] adj, boolean[] forbidden, int source) {
        // int n = g.n;
        // int[][] adj = g.adj;
        int n = adj.length;
        
        int parent[] = new int[n];
        Arrays.fill(parent, -1);
        // Create a queue for BFS
        LinkedList<Integer> queue = new LinkedList<Integer>();

        // Mark the first node as visited and enqueue it
        queue.add(source);
        parent[source] = -2;
        while (!queue.isEmpty()) {
            int s = queue.remove();

            for (int j : adj[s]) {
                if (!forbidden[j] && parent[j] == -1) {
                    parent[j] = s;
                    queue.add(j);
                }
            }
        }
        return parent;
    }

    /**
     * To sort the adjacency lists such that vertices in each list follow the given order, or increasing order when {@code permutation = null}.
     * 
     * To achive linear time, no sorting algorithm is called.
     * Instead, we build new adjacency lists by inserting their elements in the desired order.
     *
     * For example, when soritng { { 4, 1 }, { 2, 0 }, { 3, 1 }, { 2, 4 }, { 3, 0 } } with {@code permutation = null}, we have the following steps.
     * After processing 0:  { { }, { 0 }, { }, { }, { 0 } };
     * After processing 1:  { { 1 }, { 0 }, { 1 }, { }, { 0 } };
     * After processing 2:  { { 1 }, { 0, 2 }, { 1 }, { 2 }, { 0 } };
     * After processing 3:  { { 1 }, { 0, 2 }, { 1, 3 }, { 2 }, { 0, 3 } };
     * After processing 4:  { { 1, 4 }, { 0, 2 }, { 1, 3 }, { 2, 4 }, { 0, 3 } }.
     * 
     * @param adj: original adjacency lists
     * @param permutation: the specified ordering or {@code null}
     * @return the sorted adjacency lists
     **/
    public static int[][] sortAdjacencyLists(int[][] adj, int[] permutation) {
        int n = adj.length;
        if (permutation != null && permutation.length != n)
            throw new IllegalArgumentException(
                    "The order of the input graph and the size of the given permutation do not match.");
        assert(permutation == null || isPermutation(permutation));
        int[][] newlists = new int[n][];
        for (int i = 0; i < n; i++) {
            newlists[i] = new int[adj[i].length];
            // cur[i] = 0;
        }
        // System.out.println("adj = " + Arrays.deepToString(adj));

        int[] cur = new int[n]; // keep track of the position of the next neighbor of vertex i in adj[i].
        for (int i = 0; i < n; i++) {
            // for (int i = n-1; i >=0; i--) {// in the reversed order as permutation
            int v = i;
            if (permutation != null) v = permutation[i]; // the ith vertex in permutation
            for (int j : adj[v]) {
                newlists[j][cur[j]++] = v;
            }
        }
        return newlists;
    }

    /**
     * The ordering of the initial bucket for LBFS↑.
     * Sort the vertices by max{sigma^{-1}(x): x∈ N[v]}, and for vertices with the
     * same value, sort them in the reversal of their indices in sigma.
     *
     * We use a for loop to process the vertices backward according to {@code sigma}.
     * For each vertex that has not been positioned, we put it to the last available position, followed by its neighbors that have not been positioned.
     *
     * @param originalLists: the original adjacency lists, which must have been sorted with {@code sortAdjacencyLists()}.
     * @param sigma: the order of the previous search
     * @return the default order of the vertices
     **/
    private static int[] upOrdering(int[][] originalLists, int[] sigma) {
        int n = originalLists.length;
        if (sigma.length != n)
            throw new IllegalArgumentException(
                    "The order of the input graph and the size of the given permutation do not match.");
        assert(isPermutation(sigma));

        int[] newOrder = new int[n];
        boolean[] processed = new boolean[n];
        // newOrder[0] = previousOrder[n - 1];
        // processed[previousOrder[n - 1]] = true;

        int write = 0;
        for (int i = n - 1; i >= 0; i--) {
            int v = sigma[i];
            if (!processed[v]) {
                processed[v] = true;
                newOrder[write++] = v;
            }
            for (int j : originalLists[v]) { // we go through N(v) even if v itself has been positioned.
                if (!processed[j]) {
                    processed[j] = true;
                    newOrder[write++] = j;
                }
            }
        }
        return newOrder;
    }

    /**
     * The wrapper for LBFS, for which there is no previous order and the indicator for LBFS↑ is {@code false}.
     *
     * @param g - the input graph
     * @return an LBFS ordering starting from an arbitrary vertex
     */
    public static int[] LBFS(Graph g) {
        return genericLBFS(g, null, false);
    }

    /**
     * The wrapper for LBFS+, for which the indicator for LBFS↑ is {@code false}.
     *
     * @param g - the input graph
     * @param permutation - the previous order
     * @return <em>the</em> LBFS+ ordering after {@code permutation}
     */
    public static int[] LBFSplus(Graph g, int[] permutation) {
        return genericLBFS(g, permutation, false);
    }

    /**
     * The wrapper for LBFS↑.
     *
     * @param g - the input graph
     * @param permutation - the previous order
     * @return a LBFS↑ ordering after {@code permutation}
     */
    public static int[] LBFSup(Graph g, int[] permutation) {
        return genericLBFS(g, permutation, true);
    }

    /**
     * For the recognition of unit interval graphs.
     * 
     * @param g: graph
     * @param endvertices: an end vertex for each component     
     * @return an LBFS ordering that always choose a minimum-degree vertex when there are choices.
   */
    public static int[] LBFSdelta(Graph g, int[] endVertices) {
        int[] degree = Arrays.stream(g.adj).mapToInt(a -> a.length).toArray();
        for (int i : endVertices)
            degree[i] = 0; // force the end vertices to be at the beginning.
        // this permutation is sorted vertices by degree from small to large
        int[] permutation = countSortIndex(degree);
        return genericLBFS(g, permutation, false);
    }

    /**
     * A linear-time implementation of LBFS and its variants (step 2.2 of [C21]).
     *
     * The main trick is to use partition refinement instead of recording and comparing labels, which are very unlikely (if possible at all) to be done in linear time.
     * Vertices with the same lable are stored in a bucket, and buckets are organized as a linked list in decreasing order of the labels.
     * After visiting a new vertex, we split each bucket into two (if both nonempty).
     * The elements in each new bucket preserve the original related order (without sorting again).
     *
     * For example, for a five-cycle with vertex 0-1-2-3-4-
     * initially one bucket: [0, 1, 2, 3, 4] (empty label)
     * after visiting vertex 0: [1, 4] (label "0") -> [2, 3] (empty label)
     * after visiting vertex 1: [4] (label "0") -> [2] (label "1") -> [3] (empty
     * label)
     * after visiting vertex 4: [2] (label "1") -> [3] (label "4")
     * after visiting vertex 2: [3] (label "2, 4")
     *
     * @param adj           - the adjacency list of the graph
     * @param previousOrder - the previous search ordering; for LBFS↑ and LBFS+
     * @param up            - whether this is LBFS↑ search
     * @return the new LBFS order
     */
    private static int[] genericLBFS(Graph g, int[] previousOrder, boolean up) {
        int n = g.n;
        int[] newOrder = new int[n];
        int[] tieBreakers = null;
        int[] earlierNeighbors = null; // # neighbors of i that are before i in previousOrder and have not been visited in the present search. Only for LBFS↑
        int[] inversePreviousOrder = null; // For the calculation of earlierNeighbors. Only for LBFS↑

        // preprocessing for LBFS↑.
        if (up) {
            inversePreviousOrder = inversePermutation(previousOrder);
            earlierNeighbors = new int[n];
            // Arrays.fill(earlierNeighbors, 0); // Can be omitted in Java.
            for (int i = 0; i < n; i++)
                for (int j : g.adj[i])
                    if (inversePreviousOrder[j] < inversePreviousOrder[i])
                        earlierNeighbors[i]++;
        }

        /*
         * To find the next vertex of a bucket S in O(|S|) time (step 2.2 of LBFS),
         * we need the vertices in S to be sorted in a certain order.
         * More importantly, we need to maintain the order when creating a new bucket.
         * For this purpose, we preprosess the adjacency lists such that they all
         * observe the desired order.
         * For step 2.4.
         */
        int[][] adj = g.adj;
        if (previousOrder != null) {
            tieBreakers = Arrays.copyOf(previousOrder, n);
            reverse(tieBreakers); // reverse() is in-place, so we apply it to a copy of previousOrder.
            adj = sortAdjacencyLists(adj, tieBreakers); // each adjacency list follow the reversed order of previousOrder.
        }
        if (up) {
            tieBreakers = upOrdering(adj, previousOrder);
            adj = sortAdjacencyLists(adj, tieBreakers);
        }

        /*
         * It is very unlikely to implement partition refinement using existing data
         * structures in Java library.
         *
         * Each bucket, a doubly linked list, maintains the vertices with the same
         * label; and the buckets are organized in lexicographic order.
         * 
         * Initially, there is a single bucket {@code initialBucket} containing all the vertices, each in a distinct node.
         * It's 
         *   0 --> 1 --> 2 --> ...
         * if previousOrder = null; or
         *   4 --> 6 --> 3 --> 5 --> 0 --> 2 --> 1
         * if the previousOder is [1, 2, 0, 5, 3, 6, 4].
         * 
         * We'll maintain all the nodes throughout, and vertexToNode[i] refers to node containing vertex i.
         */
        Node<DoublyLinkedList<Integer>>[] vertexToBucket = new Node[n]; // Which bucket contains this vertex
        Node<Integer>[] vertexToNode = new Node[n]; // The node for the vertex
        DoublyLinkedList<Integer> initialBucket = new DoublyLinkedList<Integer>(-1);
        Node<DoublyLinkedList<Integer>> initialNode = new Node<DoublyLinkedList<Integer>>(initialBucket);
        DoublyLinkedList<DoublyLinkedList<Integer>> partition = new DoublyLinkedList<DoublyLinkedList<Integer>>(-1); // time represents when the linked list is created
        partition.insertAtBeginning(initialNode);
        for (int i = 0; i < n; i++) {
            // insertion by the previousOrder
            // previousOrder[number] = vertex
            int v = (previousOrder == null) ? i : tieBreakers[i]; // the ith of previousOrder from the last.
            Node<Integer> node = new Node<Integer>(v);
            vertexToNode[v] = node;// the node containing v, crucial for efficient implementation of step 2.4; see
                                   // comment below
            vertexToBucket[v] = initialNode; // the linked list containing the node
            // For LBFS+ and LBFS↑, elements in the bucket should be in the reversed order
            // as the previous search.
            initialBucket.insertAtEnd(node);
        }
        // initialBucket.display();

        for (int i = 0; i < n; i++) {
            /*
             * Step 2.1
             * S ← the first bucket in partition, which contains unvisited vertices with the
             * lexicographically largest label.
             */
            DoublyLinkedList<Integer> firstBucket = partition.head.element;

            Node<Integer> nodeV = null; // the vertex to be visited next
            if (up != true) {
                /*
                 * Step 2.2. v ← the last vertex of σ|S;
                 * For LBFS, any vertex in the bucket is good;
                 * For LBFSdelta, we have prepared the bucket such that the degrees are nondecreasing.
                 */
                nodeV = firstBucket.head;
                // if (DEBUG) System.out.println("i = " + i + ", v: " + nodeV.element + ": " +
                // partition);
                assert (!firstBucket.isEmpty());
            } else {
                /**
                 * This is the most nontrivial step of LBFS↑.
                 * 
                 * It needs to choose an exposed vertex if one exists and g is an interval graph
                 * (the following may make a wrong choice when g is not an interval graph).
                 * Otherwise, it chooses the vertex in the bucket that appears last in the previous search, same as LBFS+.
                 * See the paper for correctness.
                 */
                Node<Integer> nodeP = firstBucket.head, temp = firstBucket.head;
                while (temp.next != null) {
                    temp = temp.next;
                    if (inversePreviousOrder[nodeP.element] > inversePreviousOrder[temp.element])
                        nodeP = temp;
                }
                int vp = nodeP.element; // the vertex in this bucket that is the first in the previous order.

                if (earlierNeighbors[vp] > 0) {
                    nodeV = nodeP;
                } else
                    nodeV = firstBucket.head;  // The vertex in this bucket that is the last in the previous order. Recall that vertices in the bucket always keep the original order, which is reverse to the previous order.
            }

            /**
             * Step 2.3: assign i to v and remove v from the partition.
             */
            int v = nodeV.element;
            newOrder[i] = v;
            firstBucket.delete(nodeV);
            if (firstBucket.isEmpty())
                partition.delete(vertexToBucket[v]);
            vertexToNode[v] = null; // the same effect as visited = true in BFS.

            /**
             * Step 2.4: update the labels of unvisited neighors of v; adding "i".
             * For implementation, we split each bucket into two (!) buckets, neighbors of v, and non-neighbors.
             * The definition of lexicographic order ensures that new buckets from different old buckets won't mix up.
             * For example, if bucket B1 was higer than bucket B2, then B1 \cap N(v) and B1 \ N(v) are higher than B2 \cap N(v) and B2 \ N(v).
             * ! Do not create empty bucket.
             *
             * Since this step needs to be done in O(d(v)) time, we are not allowed to scan the buckets.
             * For each neighbor j of v, we use {@code vertexToBucket} and
             * {@code vertexToNode} to locate the bucket containing j and the position of j inside this bucket.
             * We remove j from this bucket and add it to a higher bucket (create one if the immediately above is not new).
             *
             * The order of vertices remaining in the bucket does not change.
             * To make sure the vertices in the new bucket are in the desired order, we insert them in the previous order (recall that we have preprocessed the adjacency list of v).
             *
             * We also need to keep track of the newly created buckets, and
             * to make sure empty buckets removed (meaning B \subseteq N[v]).
             *
             * Each j takes constant time, and the for loop takes O(d(v)) time in total.
             */
            for (int j : adj[v]) {
                Node<Integer> nodeJ = vertexToNode[j];
                if (nodeJ == null) // ignoring neighbors of v previously visited.
                    continue;
                Node<DoublyLinkedList<Integer>> nodeBucketJ = vertexToBucket[j];
                DoublyLinkedList<Integer> bucketJ = nodeBucketJ.element;

                // whether bucket list has been split when visiting v
                boolean split = (nodeBucketJ.previous != null && nodeBucketJ.previous.element.time == i);

                // we do nothing if the bucket has been not split and contains only vertex j
                if (!split && nodeJ.next == null && nodeJ.previous == null)
                    continue;

                // if |B| > 1 and B has not been split, create a new bucket and put it immediately before the current one.
                if (split == false) {
                    DoublyLinkedList<Integer> newlist = new DoublyLinkedList<Integer>(i);
                    Node<DoublyLinkedList<Integer>> newNode = new Node<DoublyLinkedList<Integer>>(newlist);
                    partition.insertBefore(nodeBucketJ, newNode);
                }

                // promote vertex j from current bucket to the newly created.
                bucketJ.delete(nodeJ);
                nodeBucketJ.previous.element.insertAtEnd(nodeJ); // the same relative order as in adj[v].
                vertexToBucket[j] = nodeBucketJ.previous; // update the bucket reference.
                // If {@code bucketJ} becomes empty, the previous must be newly created; all vertices in {@code bucketJ} are adjacent to v, and hence we set the time of the previous bucket to be i - 1.
                // Otherwise, the next bucket get the wrong value of {@code split}.  This is "lying," but the timestamp of a bucket has no other use.
                if (bucketJ.isEmpty()) {// split must be true in this case
                    nodeBucketJ.previous.element.time = i - 1;
                    partition.delete(nodeBucketJ);
                }

                // For LBFS↑, to update earlierNeighbors.
                if (up && inversePreviousOrder[j] > inversePreviousOrder[v]) {
                    earlierNeighbors[j]--;
                }
            }
        }

        return newOrder;
    }

    /**
     * The algorithm maximum cardinality search in linear time.
     *
     * @param g the input graph
     * @return an mcs ordering of g.
     **/
    public static int[] mcs(Graph g) {
        int n = g.n;
        int[] newOrder = new int[n];

        /*
         * It is very unlikely to implement partition refinement using existing data
         * structures in Java library.
         *
         * Each bucket, a doubly linked list, maintains the vertices with the same
         * label; and the buckets are organized in lexicographic order.
         * 
         * Initially, there is a single bucket {@code initialBucket} containing all the
         * vertices, each in a distinct node.
         * We'll maintain all the nodes throughout, and vertexToNode[i] refers to node
         * containing vertex i.
         * After visiting a new vertex, we split each bucket into two (if both
         * nonempty).
         */
        Node<DoublyLinkedList<Integer>>[] vertexToBucket = new Node[n]; // Which bucket contains this vertex
        Node<Integer>[] vertexToNode = new Node[n]; // The node for the vertex
        DoublyLinkedList<Integer> initialBucket = new DoublyLinkedList<Integer>(-1);
        Node<DoublyLinkedList<Integer>> initialNode = new Node<DoublyLinkedList<Integer>>(initialBucket);
        DoublyLinkedList<DoublyLinkedList<Integer>> partition = new DoublyLinkedList<DoublyLinkedList<Integer>>(-1); // time represents when the linked list is created
        partition.insertAtBeginning(initialNode);
        for (int v = 0; v < n; v++) {
            // insertion by the previousOrder
            // previousOrder[number] = vertex
            Node<Integer> node = new Node<Integer>(v);
            vertexToNode[v] = node;// the node containing v, crucial for efficient implementation of step 2.4; see
                                   // comment below
            vertexToBucket[v] = initialNode; // the linked list containing the node
            // For LBFS+ and LBFS↑, elements in the bucket should be in the reversed order
            // as the previous search.
            initialBucket.insertAtEnd(node);
        }
        // initialBucket.display();

        for (int i = 0; i < n; i++) {
            /*
             * Step 2.1
             * S ← the first bucket in partition, which contains unvisited vertices with the
             * lexicographically largest label.
             */
            DoublyLinkedList<Integer> firstBucket = partition.head.element;
            assert (!firstBucket.isEmpty());
            Node<Integer> nodeV = null; // the vertex to be visited next
            nodeV = firstBucket.head;

            /**
             * Step 2.3: assign i to v and remove v from the partition.
             */
            int v = nodeV.element;
            newOrder[i] = v;
            firstBucket.delete(nodeV);
            if (firstBucket.isEmpty())
                partition.delete(vertexToBucket[v]);
            vertexToNode[v] = null;

            /**
             * Step 2.4: update the labels of unvisited neighors of v; adding "i".
             * For implementation, we split each bucket into two (!) buckets, neighbors of
             * v, and non-neighbors.
             * The definition of lexicographic order ensures that new buckets from different
             * old buckets won't mix up.
             * For example, if bucket b1 was higer than bucket b2, then b11 and b12 (out of
             * b1) are higher than b21 and b22 (out of b2).
             * ! Do not create empty bucket.
             *
             * To achieve linear-time, this step needs to be done in O(d(v)) time, and hence
             * we are not allowed to scan the buckets.
             * For each neighbor j of v, we use {@code vertexToBucket} and
             * {@code vertexToNode} to locate the bucket containing j and the position of j
             * inside this bucket.
             * We remove j from this bucket and add it to a higher bucket.
             *
             * The order of vertices remaining in the bucket does not change.
             * To make sure the vertices in the new bucket are in the desired order, we have
             * preprocessed the adjacency list of v.
             *
             * We also need to keep track of the newly created buckets, and
             * to make sure empty buckets removed (when all the vertices in a bucket
             * promoted).
             */
            for (int j : g.adj[v]) {
                Node<Integer> nodeJ = vertexToNode[j];
                if (nodeJ == null) // ignoring neighbors of v previously visited.
                    continue;
                Node<DoublyLinkedList<Integer>> nodeBucketJ = vertexToBucket[j];
                DoublyLinkedList<Integer> bucketJ = nodeBucketJ.element;

                // whether the difference between the previous bucket and this bucket is one.
                boolean hasOnePlus = (nodeBucketJ.previous != null && nodeBucketJ.previous.element.time - nodeBucketJ.element.time == 1);

                // if the bucket contains only vertex j and the previous has a label greater than the current + 2, we increase the label
                if (!hasOnePlus && nodeJ.next == null && nodeJ.previous == null) {
                    nodeBucketJ.element.time++;
                    continue;
                }

                // otherwise, create a new bucket and put it immediately before the current one.
                if (hasOnePlus == false) {
                    DoublyLinkedList<Integer> newlist = new DoublyLinkedList<Integer>(nodeBucketJ.element.time+1);
                    Node<DoublyLinkedList<Integer>> newNode = new Node<DoublyLinkedList<Integer>>(newlist);
                    partition.insertBefore(nodeBucketJ, newNode);
                }

                // promote vertex j from current bucket to the newly created.
                bucketJ.delete(nodeJ);
                nodeBucketJ.previous.element.insertAtEnd(nodeJ); // the same relative order as in adj[v].
                vertexToBucket[j] = nodeBucketJ.previous;
                // if {@code bucketJ} becomes empty, the previous must be newly created; all
                // vertices in {@code bucketJ} are adjacent to v, and hence we set the time of
                // the previous bucket to be i - 1.
                // otherwise, the next bucket get the wrong value of split.
                if (bucketJ.isEmpty()) {// split must be true in this case
                    partition.delete(nodeBucketJ);
                }
            }
        }

        return newOrder;
    }
}
