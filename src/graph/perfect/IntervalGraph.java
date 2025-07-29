package graph.perfect;

import java.security.SecureRandom;
import java.io.*;
import java.util.*;
import java.util.stream.IntStream;

import graph.Graph;
import graph.algorithms.GraphSearch;
import static util.Permutation.*;

/**
 * 
 * @author Yixin Cao and Chenxi Liu (May 2025)
 *
 * Interval graphs

 * By default, an ordering σ is from [0..n-1] → V(G); hence, σ[i] = v meaning that vertex v is the ith visited (counted from 0).
 * See {@class graph.algorithms.GraphSearch} for more details.
 *
 */
public class IntervalGraph extends ChordalGraph {
    public static boolean DEBUG = true;// false;
    protected int[] iOrder; // an interval ordering

    /**
     * Default constructor: to generate an empty interval graph.
     * 
     * @param n: the order of the graph
     * @return an empty interval graph with n vertices.
     */
    public IntervalGraph(int n) {
        super(n);
        iOrder = IntStream.range(0, n).toArray();
    }

    /**
     * Create a random interval graph.
     * Assume that the intervals for vertices 0 to n-1 have non-decreasing left end.
     * Randomly assign the number of late neighbors for each vertex.
     * 
     * @param n: the order of the graph
     * @return a random interval graph with n vertices.
    */
    public static IntervalGraph getIntervalGraph(int n) {
        int lateNeighbors[] = new int[n];
        SecureRandom rand = new SecureRandom();
        for (int i = 0; i < n; i++) 
            lateNeighbors[i] = rand.nextInt(n - i);  // vertex i has lateNeighbors[i] neighbors after i.
        return getIntervalGraph(lateNeighbors);
    }
    
    /**
     * Generate an interval graph with the given last neighbors.
     * 
     * Add edges between i and i+1...i+lateNeighbors[i]
     * 
     * @param lateNeighbors: an int array where lateNeighbors[i] is the number of neighbors of vertex i after i.
     * @return an interval graph with n vertices.
     */
    public static IntervalGraph getIntervalGraph(int[] lateNeighbors) {
        int n = lateNeighbors.length;
        IntervalGraph g = new IntervalGraph(n);
        int[] earlyNeighbors = new int[n];
        int[] cursors = new int[n];
        for (int i = 0; i < n; i++) {
            assert(lateNeighbors[i] >= 0);
            // if (lateNeighbors[i] > n) lateNeighbors[i] = n;
            g.adj[i] = new int[lateNeighbors[i] + earlyNeighbors[i]];
            for (int j = 0; j < lateNeighbors[i]; j++) {
                g.adj[i][j] = i + j + 1;
                earlyNeighbors[i + j + 1]++;
            }
            cursors[i] = lateNeighbors[i];
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < lateNeighbors[i]; j++)
                g.adj[i + j + 1][cursors[i + j + 1]++] = i;
        }
        return g;
    }

    /**
     * An interval ordering of the graph.
     * After the first calculation, store it.
     *
     * This method assumes that {@code this} is a valid interval graph and does not check the validity of the ordering produced.
     *
     * @return an interval ordering of the graph.
     */
    public int[] getIntervalOrder() {
        if (iOrder != null)
            return iOrder;

        int[] tau = GraphSearch.LBFS(this);
        int[] tauPlus = GraphSearch.LBFSplus(this, tau);
        int[] pi = GraphSearch.LBFSup(this, tauPlus);
        iOrder = GraphSearch.LBFSplus(this, pi);
        return iOrder;
    }

    /**
     * Draw an interval model for the graph.
     *
     * TODO: use iOrder.
     *
     * TODO: output the model in tikz format.
     */
    public void visualize() {
        int clique = 0;
        for (int i = 0; i < n; i++) {
            int earlyNeighbors = 0;
            for(int j: adj[i])
                if (j < i)
                    earlyNeighbors++;
            if (earlyNeighbors+1 > clique) clique = earlyNeighbors+1;
        }
        boolean[] forbiddenColors = new boolean[clique];
        int[] colors = new int[n];
        for (int i = 0; i < n; i++) {
            for(int j: adj[i])
                if (j < i)
                    forbiddenColors[colors[j]] = true;
            int c = 0;
            while (forbiddenColors[c])
                c++;
            colors[i] = c;
            for(int j: adj[i])
                if (j < i)
                    forbiddenColors[colors[j]] = false;
        }
        if (DEBUG)
            System.out.println("coloring: " + Arrays.toString(colors));
        int w = String.valueOf(n).length();

        char[][] output = new char[clique][ (w+2)*n];
        for (int i = 0; i < n; i++) {
            int start = i * (w+2);
            output[colors[i]][start++] = '\u2517';
            String s = String.format("%0" + w + "d", i);
            for (int j = 0; j < s.length(); j++) {
                output[colors[i]][start++] = s.charAt(j);
            }
            int lastN = i;
            for(int j: adj[i])
                if (j > lastN)
                    lastN = j;
            while (start <= lastN * (w+2))
                output[colors[i]][start++] = '-';              
            output[colors[i]][start] = '\u251B';
        }
        for (char[] row: output) 
            System.out.println(new String(row));
    }

    /**
     * A linear-time algorithm to check whether <0, 1, \ldots, n-1> is an interval ordering.
     * 
     * Check for all i = 1, ..., n-1, the ith list starts from [f(i), f(i)-1, ...,
     * i+1], where f(i) is the first number in the list.
     *
     * A single triple of i, j, t with i<j<t such that
     *     sigma[i] is adjacency to sigma[t] but not sigma[j]
     * is sufficiently to refute the order.
     * 
     * For certifying algorithm, the algorithm needs to find a triple with the minimum t when sigma is not an interval ordering of g.
     * In this case, j must be t-1 (assuming g is chordal and sigma is an LBFS ordering).
     *
     * @param g: a graph (the neighborhood of each vertex must have been sorted)
     * @param badPair: if not null, the algorithm will put i and t into this array when sigma is not an interval ordering of g.
     * @return  {@code true} if <0, 1, \ldots, n-1> is an interval ordering of {@code g}.
     */
    public static boolean isIntervalOrderHelper(Graph g, int[] badPair) {
        int n = g.n;
        int[][] adj = g.adj;
        int t = n + 1, s = -1;

        for (int v = 0; v < n - 1; v++) {
            int deg = adj[v].length;
            if (deg == 0)
                continue;
            int last = adj[v][deg-1];
            if (last <= v + 1)
                continue;
            if (last - v > deg || (adj[v][deg - last + v] != v + 1)) {
                if (badPair == null)
                    return false;  // done for recognition.
                // TODO: use binary search.
                /*
                 * The following loop finds the last violation instead of the first!
                 int j = deg - 2;
                 while (j >= 0 && j > (deg - last + i) && inversedSigma[adj[v][j]] == last + 1 + j - deg) j--;
                if (last + 2 + j - deg < t) {
                    t = last + 2 + j - deg;
                    s = i;
                }
                */
                int j = 0;
                while (adj[v][j] < v)
                    j++;
                int first = j - 1;
                while (adj[v][j] == v+j-first)
                    j++;
                if (adj[v][j] < t) {
                    t = adj[v][j];
                    s = v;
                }
            }
        }
        if (s < 0)
            return true;
        badPair[0] = s;
        badPair[1] = t;
        return false;
    }

    /**
     * Is a graph an interval graph.
     * It only asks for a yes/no answer, withnot requiring a certificate.
     *
     * @param g: a graph
     * @return {@code true} if the graph is an interval graph
     */
    public static boolean isInterval(Graph g) {
        if (g.n < 4)
            return true;
        return isInterval(g, null);
    }

    /**
     * If certificate is null, it is a simple recogtion.
     * Otherwise, it returns a negative certificate.
     *
     * The certifying recognition is described in
     *     Yixin Cao, Peng Li, and Yaokun Wu, Certifying recognition of interval graphs using graph searches.
     * 
     * If g is not an interval graph and certificate is null, certificate will contains the vertex set of a hole or a Lekkerkerker--Boland graph when the method returns.
     * 
     * @param g: a graph
     * @param certificate: if not null, {@code certificate[0]} will store
     *     - an interval ordering of g when it is an interval graph;
     *     - or the vertex set of a hole or a Lekkerkerker--Boland graph otherwise.
     * @return {@code true} if the graph is an interval graph
     */
    public static boolean isInterval(Graph g, int[][] certificate) {
        boolean certifying = (certificate != null);
        if (g.n < 4)
            return true;
        int[] tau = GraphSearch.LBFS(g);
        List<Integer> vsHole = new ArrayList<Integer>();
        boolean isChordal = isLbfsPeo(g, tau, vsHole);
        if (isChordal == false) {
            certificate[0] = vsHole.stream().mapToInt(Integer::intValue).toArray();
            return false; // only proceed when g is chordal.
        }
        
        int[] tauPlus = GraphSearch.LBFSplus(g, tau);
        int[] pi = GraphSearch.LBFSup(g, tauPlus);
        int[] piPlus = GraphSearch.LBFSplus(g, pi);

        if (DEBUG)
            System.out.println("tau = " + Arrays.toString(tau)
                    + "\n tau+: " + Arrays.toString(tauPlus)
                    + "\n pi: " + Arrays.toString(pi)
                    + "\n pi+: " + Arrays.toString(piPlus));
                    
        int[] badPair = null;
        if (certifying) {
            badPair = new int[2];
        }
        Graph g2 = new Graph(GraphSearch.sortAdjacencyLists(g.subgraph(piPlus).adj, null));

        boolean ans = isIntervalOrderHelper(g2, badPair);
        if (!certifying)
            return ans;
        if (ans) {
            certificate[0] = piPlus;
        }
        else {  //  && certifying
            int[] newP = new int[g.n];
            int[] inversedPiPlus = inversePermutation(piPlus);            
            for (int i = 0; i < g.n; i++)
                newP[i] = inversedPiPlus[pi[i]];
            if (DEBUG)
                System.out.println("badPair: " + Arrays.toString(badPair) + ", n = " + g.n);
            assert(badPair[1] < g.n);
            List<Integer> vsLB = new LinkedList<Integer>();
            findCertificate(g2, newP, badPair, vsLB);
            // vsLB.forEach( (n) -> { certificate.add(piPlus[n]); } );
            certificate[0] = vsLB.stream().mapToInt(n -> piPlus[n]).toArray();            
            //            for (int i = 0; i < certificate.size(); i++) {                certificate.add(piPlus[certificate.remove(0)]);            }
        }
        return ans;
    }

    /**
     * To find a Lekkerkerker--Boland subgraph from a small (at most 10 vertices) non-interval graph
     * TODO: make sure that the first three vertices are the AT.
     * 
     * @param a graph
     * @param the vertex set of a non-interval subgraph on at most ten vertices; the redundent vertices will be removed.
     */
    private static void removeRedundant(Graph g, List<Integer> vsLB) {
        assert(vsLB.size() <= 10); // we do not use this method for a large amout of vertices
        System.out.println("vsLB: " + vsLB);
        assert(!isInterval(g.subgraph(vsLB)));
        for (int i = 0; i < vsLB.size(); i++) {
            int v = vsLB.remove(0);
            if (!isInterval(g.subgraph(vsLB))) {
                removeRedundant(g, vsLB);
                System.out.println("one vertex removed: " + v);
                return;
            }
            vsLB.add(v);
        }
    }

    /**
     * To find the Lekkerkerker--Boland graph from the asteroidal triple
     *     {q, t-1, t}
     *
     * Either p = q (step 5) or p < q && ell ~ t (step 7).
     * 
     * @param g: the new graph with vertices renumbered
     * @param t: the smallest number such that <1, 2, ..., t> is not an interval ordering
     * @param s:  
     * @param q: the smallest number in the component of G - N[t-1] containing t
     * @param vsLB: a list to store the negative certificate
     */
    private static void findLB(Graph g, int q, int s, int t, List<Integer> vsLB) {
        if (DEBUG)
            System.out.println("finding a Lekkerkerker--Boland subgraph from the AT {q, t - 1, t}."
                               + "\nq = " + q + ", s = " + s + ", t = " + t);
        int n = g.n;
        int[][] adj = g.adj; 
        vsLB.clear(); // to remove all elements added unintentionally.
        boolean[] tNeighbors = g.neighbors(t);

        // r: the first vertex adjacent to (t-1) but not t; it exists by the definition of LBFS.
        int i = 0;
        while (tNeighbors[adj[t - 1][i]])
            i++;
        int r = adj[t - 1][i];
        /*
        for (int i = 0; i < adj[t - 1].length; i++) {
            int v = adj[t - 1][i];
            if (!tNeighbors[v]) {
                r = v;
                break;
            }
        }
        */
        if (DEBUG)
            System.out.println("adj[t] = " + Arrays.toString(adj[t]) + "\nadj[t-1] = "
                    + Arrays.toString(adj[t - 1]));

        /*
         * t_q: bfs tree of G - N(t-1) rooted at q
         * Starting from t leads to a t--q path
         *     t=u1 u2 ..... ua
         */
        boolean[] forbidden = new boolean[n];
        for (int x : adj[t - 1])
            forbidden[x] = true;
        for (int j = t + 1; j < n; j++)
            forbidden[j] = true;
        int[] t_q = GraphSearch.bfsWithForbiddenVertices(g, forbidden, q);

        boolean rNeighbors[] = g.neighbors(r);

        /*
         * Case 1: r not adjacent to q
         *
         * Let ua' be the first on the t--q path that is nonadjacent to r
         * then {t-1, r, u1, ..., ua'} is a dag
         *                          t-1
         *                           |
         *             ------------r------------
         *            |       |                    |
         * t(u1) -- u2 -- u3 --  ......  -- u(a'-1) -- ua' 
         */
        if (!rNeighbors[q]) {
            if (DEBUG)
                System.out.println("{q, t - 1, t} is an AT, r not adjacent to q.");
            vsLB.add(t);
            vsLB.add(t - 1);
            vsLB.add(r);
            int v = t;
            if (DEBUG)
                System.out.println("v = " + v + ", parent[v] = " + t_q[v] + "\n" + Arrays.toString(t_q));
            while (t_q[v] > -1 && rNeighbors[t_q[v]]) {
                v = t_q[v];
                vsLB.add(v);
            }
            // after it stops, v is adjacent to r and its next (on the path from t to q) is not
            vsLB.addFirst(t_q[v]); // one of the AT, so put at the beginning
            return;
        }

        // Case 2: r adjacent to q
        vsLB.addFirst(t);
        // vsLB.add(r);

        /*
         * t_t1: bfs tree of G - N(q) rooted at (t-1)
         * Starting from t leads to a t--(t-1) path
         *     t=w1 w2 ... wb
         */
        Arrays.fill(forbidden, false);
        for (int x : adj[q])
            forbidden[x] = true;
        int[] t_t1 = GraphSearch.bfsWithForbiddenVertices(g, forbidden, t - 1);

        // u1 = w1 = t
        // u2 ~ w2 because pi+ is a perfect elimination ordering.
        int u2 = t_q[t], u3 = t_q[u2];
        int w2 = t_t1[t], w3 = t_t1[w2];

        /*
         * The rest is different from the paper.
         * We deal with a <= 4 and b <= 4 first, which allows us to handle (i) a = 3 and (iii) a = 4 and u3 ~ w4 together.
         *
         * a = b = 3
         *               t
         *            /    \
         *          u2 -- w2
         *         /   \   /  \
         * (u3) q  --  r -- t-1 (w3)
         *
         * a = 3, b = 4
         * Here u2 w3 may or may not be adjacent; when it is not, t-1 should be removed.
         *                   t
         *                /    \
         * (u3) q -- u2 -- w2 -- w3 -- t-1(w4)
         *        |     |       |       |       |
         *        --------------r-------------
         *
         * a = 4, b = 3; symmetric to above
         *
         * a = 4, b = 4; u2 = w2
         * Here u3 w3 may or may not be adjacent; when it is, r should be removed.
         *                       t
         *                       |
         * (u4) q -- u3 -- u2 -- w3 -- t-1(w4)
         *        |     |       |      |      |
         *        -------------r------------
         *
         * a = 4, b = 4; u2 != w2
         * Here edges between (u2, u3) and (w2, w3) may or may not be present;
         * TODO: to be more specific
         *
         *                           t
         *                        /    \
         * (u4) q -- u3 -- u2 -- w2 -- w3 -- t-1(w4)
         *        |     |       |      |       |       |
         *        -----------------r----------------
         *
         * TODO: Change the insertion order so that the AT remains at the beginning after removing redundant vertices.
         */
        if ((u3 == q || t_q[u3] == q) && (w3 == t-1 || t_t1[w3] == t-1)) {
            if (DEBUG)
                System.out.println("a = " + (u3 == q? 3 : 4) + ", b = " + (w3 == t-1? 3 : 4));
            vsLB.addFirst(q);
            vsLB.addFirst(t-1);
            vsLB.add(u2);
            // TODO: (u2 != w2) should suffice
            if (!vsLB.contains(w2)) vsLB.add(w2);
            if (u3 != q) vsLB.add(u3);
            if (w3 != t-1) vsLB.add(w3);
            if (DEBUG)
                System.out.println("a <= 4, b <= 4." + vsLB);
            vsLB.add(r);
            removeRedundant(g, vsLB); // to remove redundant by brute-force.
            return;
        }

        /**
         * when a = 3 or (a = 4 and u3~w4); b > 4
         *
         *                       q (u3)
         *                    /    \
         *        --------- u2 -- r ---------
         *       |                                |
         * (u1) t  -- w2 -- w3 -- ... -- wb  
         *
         *                            q (u4)
         *                            |
         *                ------------ u3 -----------
         *               |      |       |                |
         * (u1) t -- u2 -- w2 -- w3 -- ... -- w(b-1) -- wb  
         */
        if (u3 == q || (t_q[u3] == q) && g.isAdjacent(u3, t_t1[w3])) {
            if (DEBUG)
                System.out.println(u3 == q?"a = 3.": "a = 4 and u3~w4");

            int u = (u3 == q)? u2 : u3; // a = 3 or 4
            /*
            int u = u2;
            if (u3 != q) { // a = 4
                u = u3;
                vsLB.remove(1);
            }
            */
            
            vsLB.addFirst(q);
            vsLB.add(u);
            boolean uNeighbors[] = g.neighbors(u);

            int w = w2; // t_t1[t];
            while (w != t - 1 && uNeighbors[w]) {
                vsLB.add(w);
                w = t_t1[w];
            }
            vsLB.addFirst(w);
            if (u3 == q) vsLB.add(r);
            return;
        }

        // when b = 3 or (b = 4 and w3~u4)
        if (w3 == t - 1 || (t_t1[w3] == t-1) && g.isAdjacent(w3, t_q[u3])) {
            if (DEBUG)
                System.out.println(w3 == t-1?"b = 3.": "b = 4 and u4~w3");

            int w = (w3 == t-1)? w2 : w3; // b = 3
            /*
            int w = w2;
            if (w3 != t-1) { // b = 4
                w = w3;
                vsLB.remove(1);
            }
            */
            if (w3 == t-1) vsLB.add(r);
            vsLB.addFirst(t - 1);
            vsLB.add(w);
            boolean bw[] = g.neighbors(w);

            int u = u2; // t_q[t];
            while (u != q && bw[u]) {
                vsLB.add(u);
                u = t_q[u];
            }
            vsLB.addFirst(u);
            return;
        }

        /*
         * In the rest, a > 3 and b > 3.
         * We need at most five vertices from each path: u1...u5, w1...w5
         *
         * TODO: draw the diagrams
         */
        // 
        vsLB.add(r);
        i = 2;
        do {
            if (!vsLB.contains(u2))
                vsLB.add(u2);
            u2 = t_q[u2];
            i++;
        } while (i <= 5 && u2 != q);
        if (!vsLB.contains(u2))
            vsLB.addFirst(u2);

        i = 2;
        do {
            if (!vsLB.contains(w2))
                vsLB.add(w2);
            w2 = t_t1[w2];
            i++;
        } while (i <= 5 && w2 != t - 1);
        if (!vsLB.contains(w2))
            vsLB.addFirst(w2);
        removeRedundant(g, vsLB);
        /*
        int u3 = t_q[u2], u4 = t_q[u3], w3 = t_t1[w2], w4 = t_t1[w3];
        System.out.println("u3 = " + u3 + ", parent[u3] = " + t_q[u3] + ", w3 = " + w3 + ", parent[w3] = " + t_t1[w3]);
        if (u4 == q)
            certificate.addFirst(u4);
        else {
            // certificate.remove(1); //remove r (redundant)
            certificate.add(u4);
            certificate.addFirst(t_q[u4]);
        }
        if (w4 == piPlus[t-1])
            certificate.addFirst(w4);
        else {
            certificate.add(w4);
            certificate.addFirst(t_t1[w4]);
        }
        
        if (g.isAdjacent(u3, w4)) {
            certificate.add(u2);
            certificate.add(u3);
            return;
        }
        if (g.isAdjacent(u4, w3)) {
            certificate.add(w2);
            certificate.add(w3);          
            return;
        }
        certificate.add(u3);
        certificate.add(w3);
        
        if (!bw[u3])
            certificate.add(u2);
        if (!bu[w3])
            certificate.add(w2);
        if (u2 == w2)
            certificate.add(u2);
        */
        return;
    }

    /**
     * To find the Lekkerkerker--Boland graph from the asteroidal triple
     *     {p, t, z}
     *
     * @param g: the new graph with vertices renumbered
     * @param p: the largest number such that t in Sp.
     * @param q: the smallest number in the component of G - N[t-1] containing t
     * @param s: 
     * @param t: the smallest number such that <1, 2, ..., t> is not an interval ordering
     * @param z: the other vertex in the AT
     * @param certificate: a list to store the negative certificate
     */
    private static void findLB(Graph g, int p, int q, int s, int t, int z, int ell, List<Integer> vsLB) {
        if (DEBUG)
            System.out.println("finding a Lekkerkerker--Boland subgraph from the AT {p, t, z}."
            + "\np = " + p + ", q = " + q + ", s = " + s + ", t = " + t + ", z = " + z);
        assert (t < z);
        int n = g.n;
        int[][] adj = g.adj;

        vsLB.clear();
        // Vertices t and z are always in the Lekkerkerker--Boland graph, and belong to its asteroidal triple.
        vsLB.add(t);
        vsLB.add(z);

        /*
         * tbar: the first vertex in N(q) \ N[t]
         * zbar: the first vertex in N(q) \ N[z]
         * j: ..
         *
         * zbar <= tbar since t < z (11 in the paper).
         */
        int tbar = 0, zbar = 0, j = 0;
        boolean[] tNeighbors = g.neighbors(t);
        boolean[] zNeighbors = g.neighbors(z);
        int i = 0;
        while (zNeighbors[adj[q][i]])
            i++;
        zbar = adj[q][i];
        while (tNeighbors[adj[q][i]])
            i++;
        tbar = adj[q][i];
        /*
        for (i = adj[q].length - 1; i >= 0; i--) {
            int v = adj[q][i];
            if (!bt[v])
                tbar = v;
            if (!bz[v])
                zbar = v;
        }

        ii = q;
        while (ii < t && (tNeighbors[ii] || !zNeighbors[ii]))
            ii++;
        ell = ii;

          for (i = q; i < t; i++)
            if (!tNeighbors[i] && zNeighbors[i]) { ell = i; break;}
        */
        boolean[] qNeighbors = g.neighbors(q);
        i = 0;
        while (adj[tbar][i] > tbar || qNeighbors[adj[tbar][i]]) i++;
        j = adj[tbar][i];
        /*
          for (i = 0; i < adj[tbar].length; i++) {
            int v = adj[tbar][i];
            if (v < tbar && !qNeighbors[v]) {
                j = v;
                break;
            }
        }
        */
        if (DEBUG)
            System.out.println("tbar = " + tbar + ", zbar = " + zbar + ", ell = " + ell + ", j = " + j);
        
        /*
         * Case 1: there exists a common neighbor k of t and z in Sq
         *
 * 1.1. tbar == zbar (t !~ zbar)
 * 1.1.1. i not~ j and i ~ ell
 *           j
 *           |
 *         tbar (=zbar)
 *        /   \   
 *      i  --  ell
 *    /           \   
 *  t               z
 *
 * 1.1.2. i not~ ell (i in Sq and i not~ j; tbar ~ k)
 *           j
 *           |
 *         tbar (=zbar)
 *        / |  \   
 *       i   |  ell
 *     / \  | /   \   
 *    t  --  k --  z
 *
 * 1.1.3. i ~ j;
 *                  t
 *               /    \
 *             i   --  k
 *        /   |   X   |   \
 *      j-- tbar --ell -- z
 * 
 * 1.2. tbar != zbar (t ~ zbar)
 * 1.2.1. z ~ tbar (ell is not added).
 *                 j
 *              /    \
 *         zbar -- tbar
 *          /   \   /   \
 *         t  --  k  --  z
 *
 * 1.2.2. z !~ tbar;
 *                  t
 *               /    \
 *           zbar --  k
 *        /   |   X   |   \
 *      j-- tbar --ell -- z
 *
 */
        for (int k = q; k < t - 1; k++)
            if (tNeighbors[k] && zNeighbors[k]) {
                if (DEBUG)
                    System.out.println("AT {q, t, z}, there exists a common neighbor k of t and z in Sq");
                vsLB.add(j);
                vsLB.add(tbar);
                if (!zNeighbors[tbar]) 
                    vsLB.add(ell); // Always added when tbar = zbar.
/*
                if (zbar != tbar)
                    vsLB.add(zbar); // todo: check whether zbar can be the same as other vertices.

                if (tNeighbors[zbar]) { // when t ~ zbar; note that zbar != tbar
                    vsLB.add(k); 
                    // removeRedundant(g, vsLB);
                }
                else { // t is not adjacent to zbar.
                    if (DEBUG)
                        System.out.println("t is not adjacent to zbar: " + vsLB);
                    // i < t and i \in N(t) \ N[z].
                    i = 0; 
                    while (!tNeighbors[i] || zNeighbors[i]) i++;
                    vsLB.add(i);
                    // If i \in N(ell) \ N[j], then the path [t i ell z] avoids j, and k is not necessary.
                    if (g.isAdjacent(i, j) || !g.isAdjacent(i, ell))
                        vsLB.add(k);                      
                } 
 */
                if (zbar == tbar) { // t is not adjacent to zbar.
                    if (DEBUG)
                        System.out.println("t is not adjacent to zbar: " + vsLB);
                    // i < t and i \in N(t) \ N[z].
                    i = 0; 
                    while (!tNeighbors[i] || zNeighbors[i]) i++;
                    vsLB.add(i);
                    // If i \in N(ell) \ N[j], then the path [t i ell z] avoids j, and k is not necessary.
                    if (g.isAdjacent(i, j) || !g.isAdjacent(i, ell))
                        vsLB.add(k);                      
                    /*
                    for (i = 0; i < t; i++) 
                        if (tNeighbors[i] && !zNeighbors[i]) {
                            vsLB.add(i);
                            if (g.isAdjacent(i, j) || !g.isAdjacent(i, ell))
                                vsLB.add(k);
                        }
                    */
                }
                else  { // when t ~ zbar; note that zbar != tbar
                    vsLB.add(zbar);
                    vsLB.add(k); 
                    // removeRedundant(g, vsLB);
                }
                return;                
            }

        boolean forbidden[] = new boolean[n];
        Arrays.fill(forbidden, true);
        for (i = q; i < t; i++) {
            forbidden[i] = false;
        }
        forbidden[z] = false;
        int t_z[] = GraphSearch.bfsWithForbiddenVertices(g, forbidden, t);
        if (DEBUG)
            System.out.println("before Case 2.1: parent: " + Arrays.toString(t_z) + ", forbidden: " + Arrays.toString(forbidden));
        
        /*
         * Case 2. There is a t--z path using vertices in Sq.
 *         
 * 2.1. z !~ tbar
 *            j
 *            |
 *          tbar
 *          /   \
 *   t --   ...   --  z
 *            Sq
 *
 * 2.2. z ~ tbar (tbar != zbar).
 *                 j
 *              /    \
 *        zbar  --  tbar
 *        /   X     X   \
 *       t --   ...   --  z
 *               Sq
 *
 */
        if (t_z[z] >= 0) {
            if (DEBUG)
                System.out.println("minimize2, case 2.1");
            vsLB.add(j);
            vsLB.add(tbar);
            int v = t_z[z];
            while (v != t) {
                vsLB.add(v);
                v = t_z[v];
            }
            if (zNeighbors[tbar]) // z ~ tbar
                vsLB.add(zbar);
            return;
        }

        /*
         * Case 3. There is no t--z path using vertices in Sq.
 * 3.1. z !~ tbar
 *              u3
 *               |
 *             tbar
 *          /   |   \
 *   t -- s -- v -- ell -- z
 *
 *              u4
 *               |
 *              u3
 *               |
 *  t -- s -- u2 -- ell -- z     (u2 can be replaced by tbar?)
 *
 * 3.2. z ~ tbar (tbar != zbar).
 *                u3
 *             /      \
 *        zbar  --  tbar
 *       /   X   X   X    \
 *     t -- s -- v -- ell -- z
 *
 *              u4
 *               |
 *              u3
 *           /      \
 *        zbar -- tbar
 *         /           \
 *        t             z
 */
        for (i = 0; i < t; i++) 
            forbidden[i] = false;
        for (i = t; i < n; i++) 
            forbidden[i] = true;
        for (int x : adj[t])
            forbidden[x] = true;
        t_z = GraphSearch.bfsWithForbiddenVertices(g, forbidden, p);
        if (DEBUG)
            System.out.println("parent: " + Arrays.toString(t_z) + ", forbidden: " + Arrays.toString(forbidden));
        
        int u2 = t_z[t - 1];
        int u3 = t_z[u2];
        vsLB.add(u3);
        vsLB.add(s);
        vsLB.add(tbar);
        if (zbar != tbar)
            vsLB.add(zbar);
        vsLB.add(ell); 

        boolean[] b = g.neighbors(u3);
        for (i = 0; i < adj[q].length; i++) {
            int v = adj[q][i];
            if (v < q && !b[v]) {
                vsLB.add(v);
                removeRedundant(g, vsLB);
                return;
            }
        }
        vsLB.addFirst(t_z[u3]);
        removeRedundant(g, vsLB);        
    }

    /**
     * TODO: return an asteroidal triple (not necessarily for a Lekkerkerker--Boland graph) when vsLB==null.
     */
    private static void findCertificate(Graph g, int[] pi, int[] badPair, List<Integer> vsLB) {
        int n = g.n;
        int[][] adj = g.adj;
        int[] inversedPi = inversePermutation(pi);
        int s = badPair[0];
        int t = badPair[1];
        int[] at = new int[3];

        boolean forbidden[] = new boolean[n];
        for (int x : adj[t])
            forbidden[x] = true;
        for (int j = t; j < n; j++)
            forbidden[j] = true;
        int parent1[] = GraphSearch.bfsWithForbiddenVertices(g, forbidden, t - 1);

        /*
         * p is the last such that p < t and t in Sp
         * We search for it using Lemma 3.2 instead of the definition.
         */
        int i = 0;
        while (i < n && parent1[i] == -1)
            i++;
        int p = i;

        for (int x : adj[t])
            if (x < t)
                forbidden[x] = false;
        for (int x : adj[t - 1])
            forbidden[x] = true;
        int parent2[] = GraphSearch.bfsWithForbiddenVertices(g, forbidden, t);
        i = 0;
        while (i < n && parent2[i] == -1)
            i++;
        int q = i;
        // System.out.println("p = " + p + ", q = " + q);
        // int p = minimumReach(g, sigma, t - 1, t);
        // int q = minimumReach(g, sigma, t, t - 1);
        if (p == q) {
            findLB(g, q, s, t, vsLB);
            return;
        }

        /*
         *                               |-------Sq--------|
         * 0, ... [ p, p+1, ..., [ q, q+1, ..., t-1 ], t, ... ], .... n-1
         *          |--------------------Sp---------------------|         
         */
        int left = q, right = q;  
        for (i = q + 1; i < t; i++)
            if (inversedPi[i] > inversedPi[right])
                right = i;
        boolean[] f = new boolean[n];
        for(int j: adj[right])
            f[j] = true;
        for (i = q + 1; i < t; i++)
            if (!f[i] && inversedPi[i] < inversedPi[left])
                left = i;
        int z = n;
        for (i = 0; i < adj[left].length; i++) {
            int v = adj[left][i];
            if (v >= t && v < z)
                z = v;
        }
        if (z == t) {
            findLB(g, q, s, t, vsLB);
        } else {
            /*
            // find ell: the vertex in Sq\N[q] that appears the first in pi.
            int ell = q;
            boolean[] qNeighbors = new boolean[n];
            for (int x : adj[q])
                qNeighbors[x] = true;
            for (i = q + 1; i < t; i++) {
                if (!qNeighbors[i] && inversedPi[i] < inversedPi[ell])
                    ell = i;
            }
            */
            findLB(g, p, q, s, t, z, left, vsLB);
        }
    }

    public static void main(String[] s) {
        if (s.length < 2 || (s[0].charAt(1) != 'g' && s[0].charAt(1) != 't')) {
            System.out.println("Usage:"
                               + "\n    java graph.perfect.IntervalGraph -g <file-name>, to generate a random interval graph and write it to <file-name>;"
                               +"\n    java graph.perfect.IntervalGraph -t <file-name>, to test whether the graph stored in <file-name> is an interval graph.");
            return;
        }
        
        if (s[0].charAt(1) == 'g') {
            IntervalGraph.getIntervalGraph(10).writeToDIMACS(s[1]);
            return;
        }
        
        // List<Integer> certificate = new LinkedList<Integer>();
        Graph g = new Graph(s[1]);
        System.out.println("The input graph: " + g);
        int[][] certificate = new int[1][];
        if (isInterval(g, certificate)) {
            IntervalGraph ig = new IntervalGraph(g.n);
            ig.adj = g.adj;
            ig.iOrder = certificate[0];
            System.out.println("The graph is an interval graph.  The following interval representation uses unicode characters, and hence is not readable if your system does not support them.");
            ig.visualize();
        }
        else
            System.out.println("The graph is not an interval graph, the vertex set of the found Lekkerkerker--Boland graph: " + Arrays.toString(certificate[0]));
    }
}
