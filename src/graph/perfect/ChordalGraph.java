package graph.perfect;

import java.util.*;
import graph.Graph;
import graph.algorithms.GraphSearch;
import static util.Permutation.*;

/**
 * 
 * @author Yixin Cao and Chenxi Liu (May 2025)
 *
 * Chordal graphs
 * 
 */
public class ChordalGraph extends Graph {
    public static boolean DEBUG = false;// true;
    protected int[] pOrder; // a perfect elimination ordering
    
    public ChordalGraph(int n) {
        super(n);
    }
    
    /**
     * Construct an Erdős–Rényi graph with probability p, and then do random fill in.
     *
     * The probability should be set significantly smaller than desired, because there would be a large number of fill-in edges.
     *
     * TODO: use a randome PEO.
     */
    public ChordalGraph(int n, double p) {
        super(n);
        boolean[][] matrix = new boolean[n][n];
        //for (int i = 0; i < n; i++) 
        //for (int j = i + 1; j < n; j++) 
        //     if (p >= Math.random()) matrix[i][j] = true;
 
       for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) 
                if (p >= Math.random()) matrix[i][j] = true;
            for (int j = i + 1; j < n; j++) {
                for (int k = j + 1; k < n; k++)
                    if (matrix[i][j] && matrix[i][k])
                        matrix[j][k] = matrix[k][j] = true;
            }
        }
        // shuffle the vertex set
        int[] permutation = shuffle(n);
        boolean[][] shuffledMatrix = new boolean[n][n];
        for (int i = 0; i < n; i++) 
            for (int j = 0; j < n; j++) 
                shuffledMatrix[permutation[i]][permutation[j]] = matrix[i][j];
        this.adj=fromAdjacencyMatrix(shuffledMatrix);
    }
    
    /**
     * decide whether the graph is chordal.
     */
    public static boolean recognize(Graph graph) {
        return false;
    }

    /**
     * find a maximum clique
     */
   public int[] maxClique() {
        return null;
    }

    /**
     * Calculate the chromatic number of the graph.
     */
   public int chromatic() {
       return maxClique().length;
   }
    
    /**
     * Find a minimum coloring of the graph.
     * ans[i] is the set of vertices receiving color i.
     */
   public int[][] coloring() {
       
       return null;
    }

    /**
     * Is a graph chordal.
     */
    public static boolean isChordal(Graph g) {
        return isChordal(g, null);
    }

    /**
     * The vertex set of a hole will be stored *orderly* in certificate.
     */
    public static boolean isChordal(Graph g, List<Integer> vsHole) {
        int[] tau = GraphSearch.LBFS(g);
        if (DEBUG)
            System.out.println("graph: " + g + "LBFS: " + Arrays.toString(tau));
        return isLbfsPeo(g, tau, vsHole);
    }

    /**
     * A linear-time algorithm to check whether an LBFS ordering {@code vertexOrder} of g is a perfect elimination ordering.
     *
     * Phase 1: find the smallet number i such that G[ {sigma[0], ..., sigma[i]} ] is not chordal.
     * Use the ideas from Gustavson (1972), Whitten (1978), Tarjan&Yannakakis 1984
     *
     * Two values for each vertex,
     * f(v): the last vertex in N_σ(v) (or v when none has been met)
     * index(v): the index of the mininum vertex in N_σ[v] that has been processed.      
     * 
     * Process from n-1 to 0.
     * Initially, f[v] = v; index[sigma(i)] = i.
     * The first time we process an earlier neighbor w of v, f(v) <- w.
     * Every time we process an earlier neighbor w adjacent to v, index (v) <- a^{-1}(w).
     *
     * Phase 2: find a hole.
     * Let a, b, c be the bad triple found above, with b--a--c.
     * We search for an induced b--c path in G - N[a].
     *
     * @param g: a graph
     * @param sigma: an ordering of the vertex set of g
     * @return true if sigma is a perfect elimination ordering of g.
     */
    public static boolean isLbfsPeo(Graph g, int[] vertexOrder, List<Integer> certificate) {
        int n = g.n;
        int[][] adj = g.adj;
        // int[][] adj = GraphSearch.sortAdjacencyLists(g.adj, vertexOrder);
        int[] p = inversePermutation(vertexOrder);
        int[] badTriple = null;

        // See comments above.
        int[] f = new int[n], ind = new int[n];
        //for (int i = 0; i < n - 1; i++) {
        for (int i = n - 1; i >= 0; i--) {
            int v = vertexOrder[i];
            f[v] = v; // vertex
            ind[v] = i; // number
            for (int w: adj[v]) { // neighbor
                if (p[w] < i)
                    continue; // only care about latter neighbors.
                ind[w] = i;
                if (f[w] == w)
                    f[w] = v;
            }
            for (int w: adj[v]) { // neighbor
                if (p[w] < i || ind[f[w]] <= i)
                    continue; // only care about latter neighbors.
                if (certificate == null) 
                    return false; //only care about yes/no.
                int[] t = { w, v, f[w] };
                badTriple = t;  // find the bad triple minimizing w.
            }
        }

        if (badTriple == null)
            return true;
        boolean[] forbidden = new boolean[n];
        for (int j : adj[badTriple[0]]) 
            forbidden[j] = true;
        forbidden[badTriple[0]] = true;
        forbidden[badTriple[2]] = false;
        int parent[] = GraphSearch.bfsWithForbiddenVertices(g, forbidden, badTriple[1]);

        certificate.add(badTriple[0]);
        int v = badTriple[2];
        do {
            certificate.add(v);
           // System.out.println("v: " + v + ", parent[v]" + parent[v]);
            v = parent[v];
        } while (v >= 0);
        return false;
    }

    /**
     * Whether does a subset of vertices induce a hole, i.e., a simple cycle of length at least four. 
     * 
     * @param g: a graph
     * @param vertexSet: a subset of vertices in {@code g}
     * @return {@code true} if the subgraph induced by {@code vertexSet} is a hole.
     */
    public static boolean isHole(Graph g, int[] vertexSet) {
        return (vertexSet.length > 3 && isInducedCycle(g, vertexSet));
    }

    /**
     * Whether does a subset of vertices induce a cycle.
     *
     * For each vertex v, we check the intersection of N(v) and {@code vertexSet}.  There must be two of them, and they must be the circular predecessor and the successor of v.
     * 
     * @param g: a graph
     * @param vertexSet: a subset of vertices in {@code g}
     * @return {@code true} if the subgraph induced by {@code vertexSet} is a cycle.
     */
    public static boolean isInducedCycle(Graph g, int[] vertexSet) {
        int n = g.n, h = vertexSet.length;
        if (h < 3)
            throw new IllegalArgumentException("There must be at least three vertices.");
        /**
         * The array {@code reverse} is to store the position of each vertex on the hole.  The value {@code reverse[i]} is
         *     --  0 if vertex i is not in {@code vertexSet};
         *     --  (the position of vertex i in {@code vertexSet}) + 1 otherwise.
         * The reverse mapping will be a -> [1..|a|].
         *
         * An alternative approach is to use -1 to mark other vertices:
         * Arrays.fill(reverse, -1);
         */
        int[] reverse = new int[n];
        for (int i = 0; i < h; i++)
            reverse[vertexSet[i]] = i + 1;

        for (int i = 0; i < h; i++) {
            int v = vertexSet[i], count = 0;
            for (int j = 0; j < g.adj[v].length; j++) {
                int x = g.adj[v][j]; // neighor of v
                if (reverse[x] == 0)
                    continue; // skip x if x is not in a.
                int diff = (reverse[x] - reverse[v] + h) % h; // must be 1 or -1.
                if (diff != 1 && diff != h - 1)
                    return false;
                count++;
            }
            // There must be precisely two neighbors: otherwise, it may consists of a bunch of paths.
            if (count != 2)
                return false;
        }
        return true;
    }

    /**
     * To check whether an MCS ordering {@code vertexOrder} of g is a perfect elimination ordering.
     *
     * While MCS is simpler than LBFS, verifying MCS orderings is more complicated. 
     * Method {@code isLbfsPeo} does not work for MCS orderings, which needs a special violating triple.  
     *
     */
    public static boolean isMcsPeo(Graph g, int[] vertexOrder, List<Integer> vsHole) {
        throw new UnsupportedOperationException("Not implemented yet");
    }
    
    /**
     * Use MCS to check whether a graph is chordal.
     * 
     * @param g
     * @return
     */
    public static boolean isChordalMCS(Graph g, List<Integer> vsHole) {
        int[] order = GraphSearch.mcs(g);
        System.out.println("graph: " + g);
        System.out.println("tau: " + Arrays.toString(order));
        return isMcsPeo(g, order, vsHole);
    }

    public static boolean isChordalMCS(Graph g) {
        return isChordalMCS(g,  null);
    }

    /**
     * A dag graph of order n.
     * It is the net when n = 6;
     *
     *            0
     *            |
     *            1
     *         / | \
     *    2--3...n-2--n-1     
     */
    public static ChordalGraph net(int n) {
        ChordalGraph g = new ChordalGraph(n);
        g.adj[0] = new int[1]; g.adj[0][0] = 1;
        g.adj[2] = new int[1]; g.adj[2][0] = 3;
        g.adj[n-1] = new int[1]; g.adj[n-1][0] = n-2;
        g.adj[1] = new int[n - 3]; g.adj[1][0] = 0;
        for (int i = 3; i < n - 1; i++) {
            g.adj[1][i-2] = i;
            g.adj[i] = new int[3];
            g.adj[i][0] = 1;
            g.adj[i][1] = i - 1;
            g.adj[i][2] = i + 1;
        }        
        return g;
    }

    /**
     * A ddag graph of order n.
     * It is the sun when n = 6;
     *
     *            0
     *         /    \
     *       1  --   2
     *     /  \ X  /   \
     *    3--4...n-2--n-1     
     */
    public static ChordalGraph sun(int n) {
        ChordalGraph g = new ChordalGraph(n);
        g.adj[0] = new int[2]; g.adj[0][0] = 1; g.adj[0][1] = 2;
        g.adj[3] = new int[2]; g.adj[3][0] = 1; g.adj[3][1] = 4;
        g.adj[n-1] = new int[2]; g.adj[n-1][0] = 2; g.adj[n-1][1] = n-2;
        g.adj[1] = new int[n - 2]; g.adj[1][0] = 0; g.adj[1][1] = 2; g.adj[1][2] = 3;
        g.adj[2] = new int[n - 2]; g.adj[2][0] = 0; g.adj[2][1] = 1; g.adj[2][n - 3] = n-1;
        for (int i = 4; i < n - 1; i++) {
            g.adj[1][i-1] = g.adj[2][i-2] = i;
            
            g.adj[i] = new int[4];
            g.adj[i][0] = 1; g.adj[i][1] = 2;
            g.adj[i][2] = i - 1;
            g.adj[i][3] = i + 1;
        }        
        return g;
    }

    /**
    * Legacy code, not in use any more. To be deleted.
    */
    public static void main(String[] args) {
        if (args.length == 0) {
            System.out.println(
                    "Usage: program_name <Interval | UnitInteval |Chordal |dInterval | dUnitInteval | dChordal> <t> <n>: generate t Interval/UnitInteval/Chordal graphs ('d' indicates disconnected) of order n. \n"
                            +
                            "program_name <Interval | UnitInteval |Chordal |dInterval | dUnitInteval | dChordal> <t> <n1> <n2>: generate t Interval/UnitInteval/Chordal graphs ('d' indicates disconnected) of order in the range between n1 and n2 (n2 >= n1). \n");
            return;
        }
        /*      Scanner keyboard = new Scanner(System.in);
        String type = keyboard.next().toLowerCase();
        char c = type.charAt(0);
        boolean connected = (c != 'd');
        if (!connected) c = type.charAt(1);
        
        int t = keyboard.nextInt(); // number of graphs
        int n1 = keyboard.nextInt();
        int n2 = n1;
        if (keyboard.hasNextInt())
            n2 = keyboard.nextInt();
        */

        String type = args[0].toLowerCase();
        char c = type.charAt(0);
        boolean connected = (c != 'd');
        if (!connected) c = type.charAt(1);
        
        int t = Integer.parseInt(args[1]); // number of graphs
        int n1 = Integer.parseInt(args[2]);
        int n2 = n1;
        if (args.length > 3)
            n2 = Integer.parseInt(args[3]);
        String path = "./";
        String graphName = "chordal";
        
        for (int i = 0; i < t; i++) {
            String filename = path + graphName + i + ".txt";
            Graph g = new ChordalGraph(n1);
            if (DEBUG)
                System.out.println(Arrays.deepToString(g.adj));
            g.writeToDIMACS(filename);
            System.out.println(IntervalGraph.isInterval(g));
        }
        return;
    }
}
