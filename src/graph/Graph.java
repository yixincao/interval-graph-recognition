package graph;

import java.io.*;
import java.util.*;
// import java.util.stream.IntStream;

/**
 * @author Yixin Cao and Chenxi Liu (May 2025)
  *
 * Immutable graphs - adjacency list representation
 *
 * For simplicity, we use a two-dimensional array as the internal storage. This only works in Java.
 * 
 */
public class Graph {
    // TODO: hide the internal data
    // TODO: use a Vertex class to avoid confusion with numbers.
    public int n;
    public int[][] adj;

    /*
     *        2
     *       / \
     * 0--1--3--4
     */
    private static boolean[][] bull = {
        { false, true, false, false, false },
        { true, false, true, true, false },
        { false, true, false, true, false },
        { false, true, true, false, true },
        { false, false, false, true, false }
    };
    public static Graph BULL = new Graph(bull);

    /*
     *         0
     *         |
     *         1
     *       /  \
     * 4--2 -- 3--5
     */
    private static int[][] net = {
        { 1 },
        { 0, 2, 3 },
        { 1, 3, 4 },
        { 1, 2, 5 },
        { 2 },
        { 3 }
    };
    public static Graph NET = new Graph(net);

    /*
     *          1
     *        /  \
     *       4 -- 5
     *      /  \/  \
     *    3 -- 0 -- 2     
     */
    private static int[][] sun = {
        { 2, 3, 4, 5 },
        { 4, 5 },
        { 0, 5 },
        { 0, 4 },
        { 0, 1, 3, 5 },
        { 0, 1, 2, 4 }
    };
    public static Graph SUN = new Graph(sun);

    /*
     *     2
     *     |
     * 1--0--3
     */
     private static int[][] claw = {
        { 1, 2, 3 },
        { 0 },
        { 0 },
        { 0 }
    };
    public static Graph CLAW = new Graph(claw);

    /*
     *          6
     *          |
     *    2 -- 3 -- 4
     *    |  \ | /  |
     *    1 -- 0 -- 5     
     */
    private static int[][] top = {
        { 1, 2, 3, 4, 5 },
        { 0, 2 },
        { 0, 1, 3 },
        { 0, 2, 4, 6 },
        { 0, 3, 5 },
        { 0, 4 },
        { 3 }
    };
    public static Graph WHIPPING_TOP = new Graph(top);

    /*
     *              6
     *              |
     *              3
     *              |
     * 4 -- 1 -- 0 -- 2 -- 5
     */
     private static int[][] longClaw = {
        { 1, 2, 3 },
        { 0, 4 },
        { 0, 5 },
        { 0, 6 },
        { 1 },
        { 2 },
        { 3 }
    };
    public static Graph LONG_CLAW = new Graph(longClaw);

    private int edgeNum; // not used in the present project.
    private String[] labels; // not used in the present project.
    private String comment; // "interavl graph" etc.

    /**
     * Construct a graph out of a text file.
     */
    public Graph(String path) {
        try {
            readFile(path);
        } catch (FileNotFoundException e) {
            //            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }
    
    /**
     * Construct a graph on n vertices and 0 edges.
     */
    public Graph(int n) {
        this.n = n;
        adj = new int[n][0];
        for (int i = 0; i < n; i++) 
            adj[i] = new int[0];
    }

    /**
     * Construct an Erdős–Rényi graph with probability p.
     */
    public Graph(int n, double p) {
        this.n = n;
        boolean[][] matrix = new boolean[n][n];
        for (int i = 0; i < n; i++) 
            for (int j = i + 1; j < n; j++) 
                if (p >= Math.random()) matrix[i][j] = true;
        adj = fromAdjacencyMatrix(matrix);
    }
    
    /**
     * Construct a graph with given adjacency lists.
      */
    public Graph(int[][] adj) {
        this(adj.length);
        this.adj = adj;
    }
     
    /**
     * Construct a graph out of a boolean adjacency matrix.
     *
     * This is mainly for testing purpose: it is not linear in general.
     */
    public Graph(boolean[][] matrix) {
        this(matrix.length);
        adj = fromAdjacencyMatrix(matrix);
    }
    /*
      @SuppressWarnings("unchecked")
      public Graph(boolean[][] matrix) {
      n = matrix.length;
      adj = new int[n][];
      
      List<Integer>[] adjList = new ArrayList[n];
      for (int i = 0; i < n; i++) 
      adjList[i] = new ArrayList<Integer>();
      for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) 
      if (matrix[i][j]) {
      adjList[i].add(j);
      adjList[j].add(i);
      }
      adj[i] = adjList[i].stream().mapToInt(Integer::intValue).toArray();
      }
      }
      */
    
    /**
     * Construct a graph on n vertices and 0 edges.
     */
    public boolean validate() {
        for (int i = 0; i < n; i++) 
            adj[i] = new int[0];

        boolean[][] matrix = new boolean[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < adj[i].length; j++)
                matrix[i][adj[i][j]] = true;
        }
        
        for (int i = 0; i < n; i++) 
            for (int j = 0; j < n; j++)
                if (matrix[i][j] != matrix[j][i]) return false;
        return true;
    }
    
    @SuppressWarnings("unchecked")
    public static int[][] fromAdjacencyMatrix(boolean[][] matrix) {
        int n = matrix.length;
        int[][] adj = new int[n][];

        List<Integer>[] adjList = new ArrayList[n];
        for (int i = 0; i < n; i++) 
            adjList[i] = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) 
                if (matrix[i][j]) {
                    adjList[i].add(j);
                    adjList[j].add(i);
                }
            adj[i] = adjList[i].stream().mapToInt(Integer::intValue).toArray();
        }
        return adj;
    }

    /*
     * returen the boolean adjacency matrix of the graph
     */
    public boolean[][] getAdjacencyMatrix() {
        boolean[][] matrix = new boolean[n][n];
       for (int i = 0; i < n; i++) {
            for (int j = 0; j < adj[i].length; j++)
                matrix[i][adj[i][j]] = matrix[adj[i][j]][i] = true;
        }
        return matrix;
    }

    /**
     * the degree of vertex v.
     */
    public int degree(int v) {
        return adj[v].length;
    }

    public String toString() {
        return Arrays.deepToString(adj);
    }

    /**
     * Construct a subgraph induced by the vertexSet.
     */
    public Graph subgraph(int[] vertexSet) {
        if (vertexSet.length == 0)
            throw new IllegalArgumentException("Empty graphs are not supported.");
        int small = vertexSet.length;
        Graph subgraph = new Graph(small);
        int[] index = new int[n];
        Arrays.fill(index, -1);
        for (int i = 0; i < small; i++) 
            index[vertexSet[i]] = i;
        for (int i = 0; i < small; i++) {
            List<Integer> list = new LinkedList<Integer>();
            for (int j : adj[vertexSet[i]])
                if (index[j] > -1)
                    list.add(index[j]);
            subgraph.adj[i] = list.stream().mapToInt(Integer::intValue).toArray();
        }
        // System.out.println("vertexSet (subgraph): " + Arrays.toString(vertexSet) + "subgraph (subgraph): " + subgraph);
        return subgraph;
    }

    public boolean isAdjacent(int u, int v) {
        // if d(u) > d(v) swap them.
        return Arrays.stream(adj[u]).anyMatch(x -> x == v);
    }

    public boolean[] neighbors(int v) {
        boolean[] vNeighbors = new boolean[n];
        for (int x : adj[v])
            vNeighbors[x] = true;
        return vNeighbors;
    }
    
    /**
     * Construct a subgraph induced by the vertexSet.
     */
    public Graph subgraph(List<Integer> vertexSet) {
        return subgraph(vertexSet.stream().mapToInt(Integer::intValue).toArray());
    }
    
    // Todo:
    public LinkedList<Graph> connectedComponents;
    public int[] nodelist; // 0...n-1 maybe not consecutive, not ordered
    public int[] simplified_nodelist; // 0...vertexNum-1 consecutive

    /**
     * Load the graph from a text file in the DIMACS format.
     * 
     * The graph file format:
     * https://users.aalto.fi/~tjunttil/bliss/fileformat.html
     * 
     * Comments:
     * In the beginning of the file, there can be comment lines that start with the
     * character c.
     * For instance {c The constraint graph of a CNF formula}.
     * 
     * Problem definition line:
     * After the comment lines, there must be the “problem definition line” of the
     * form {p edge N E}.
     * Here N is the number of vertices and E is the number of edges in the graph.
     * In the file format, the vertices are numbered from 1 to N.
     * 
     * Vertex colors:
     * After the problem definition line, the next lines of the form {n v c} define
     * the colors of the vertices
     * v is the number of a vertex and c is the color of that vertex.
     * The color c should be a non-negative integer fitting in the domain of the
     * unsigned int type in the C++ programming language.
     * It is not necessary to include a color definition line for each vertex, the
     * default color for a vertex is 0.
     * If the color of a vertex is defined more than once, the last definition
     * applies.
     * 
     * Edges. Following the color definition lines, the next E lines of the form {e
     * v1 v2}
     * describe the edges in the graph, where 1 ≤ v1,v2 ≤ N are the numbers of the
     * vertices connected by the edge.
     * Multiple definitions of the same edge are ignored.
     * 
     * @param the path of the graph file
     * @return a Graph
     * @throws FileNotFoundException when the file cannot be found.
     */
    @SuppressWarnings("unchecked")
    public void readFile(String path) throws FileNotFoundException {
        File file = new File(path);
        Scanner sc = new Scanner(file);
        List<Integer>[] adjList = null;

        /*
          this.nodelist = IntStream.range(0, n).toArray();
          this.simplified_nodelist = this.nodelist;
        */
        
        while (sc.hasNextLine()) {
            String curline = sc.nextLine();
            String[] lineStrings = curline.split(" ");
            // ignore comments the vertex colors
            
            if (lineStrings[0].charAt(0) == 'p') {
                n = Integer.parseInt(lineStrings[2]);
                edgeNum = Integer.parseInt(lineStrings[3]);
                adjList = new ArrayList[n];
                for (int i = 0; i < n; i++) 
                    adjList[i] = new ArrayList<Integer>();
            }

            if (lineStrings[0].charAt(0) == 'e') {
                assert adjList != null; // there must be the "problem definition line" of the form {p edge N E}.
                int vertex1 = Integer.parseInt(lineStrings[1]);
                int vertex2 = Integer.parseInt(lineStrings[2]);
                adjList[vertex1].add(vertex2);
                adjList[vertex2].add(vertex1);
            }
        }
        sc.close();

        this.adj = new int[n][];
        for (int i = 0; i < n; i++) {
            this.adj[i] = adjList[i].stream().mapToInt(Integer::intValue).toArray();
        }

        /*
        // process connected components
        this.connectedComponents = new LinkedList<>();
        // this.connectedComponents.add(this);
        // System.out.println(this.connectedComponents.size());
        splitConnectedComponents();
        */
    }

    /**
     * It overwrites without warning.
     */
    @SuppressWarnings("unchecked")
    public void writeToDIMACS(String filename) {
        int size = 0;
        for (int i = 0; i < n; i++)
            size += adj[i].length;
        size /= 2; // m = sum of degrees / 2

        try (FileWriter writer = new FileWriter(filename)) {
            // Write header comments
            if (comment!=null)
                writer.write("c" + comment  + " \n");
            if (labels != null) {
                writer.write("c Vertex mapping:\n");
                for (int i = 0; i < n; i++) 
                    writer.write("c " + i + " = " + labels[i] + "\n");
            }
            
            // Write problem line
            writer.write("p edge " + n + " " + size + "\n");
                
            // Write edges
            for (int i = 0; i < n; i++) 
                for (int j : adj[i])
                    if (i < j) writer.write("e " + i + " " + j + "\n");
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }
}
