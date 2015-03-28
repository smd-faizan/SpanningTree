/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cytoscape.spanningtree.internal.spanningtree;

import java.util.ArrayList;
import java.util.List;
import org.cytoscape.spanningtree.internal.CyActivator;
import org.cytoscape.model.*;
import org.cytoscape.view.model.CyNetworkView;

/**
 *
 * @author smd.faizan@gmail.com
 */
public class SpanningTreeThread extends Thread {

    public CyNetwork currentnetwork;
    public CyNetworkView currentnetworkview;
    boolean isAlternative;
    boolean isMinimum;
    String edgeWeightAttribute;

    public SpanningTreeThread(CyNetwork currentnetwork, CyNetworkView currentnetworkview, boolean isAlternative, boolean isMinimum, String edgeWeightAttribute) {
        this.currentnetwork = currentnetwork;
        this.currentnetworkview = currentnetworkview;
        this.isAlternative = isAlternative;
        this.isMinimum = isMinimum;
        this.edgeWeightAttribute = edgeWeightAttribute;
    }

    // kruskals algo
    @Override
    public void run() {
        List<CyNode> nodeList = currentnetwork.getNodeList();
        int totalnodecount = nodeList.size();
        CyTable edgeTable = currentnetwork.getDefaultEdgeTable();
        CyTable nodeTable = currentnetwork.getDefaultNodeTable();

        if (!isAlternative) {
            // first spanning tree
            double[][] adjacencyMatrixOfNetwork = createAdjMatrix(currentnetwork, nodeList, edgeTable, totalnodecount, edgeWeightAttribute);
            //printMatrix(adjacencyMatrixOfNetwork, "initail matrix");
            double[][] SpTreeAdjMatrix = createSpTreeAdjMatrix(adjacencyMatrixOfNetwork, nodeList, edgeTable, totalnodecount);
            //printMatrix(SpTreeAdjMatrix, "first spanning tree");
            createNetwork(SpTreeAdjMatrix, nodeList, nodeTable, totalnodecount);
        } else {
            // alternative sp tree
            double[][] adjacencyMatrixOfNetwork = createAdjMatrix(currentnetwork, nodeList, edgeTable, totalnodecount, edgeWeightAttribute);
            //printMatrix(adjacencyMatrixOfNetwork, "initail matrix");
            double[][] SpTreeAdjMatrix = createSpTreeAdjMatrix(adjacencyMatrixOfNetwork, nodeList, edgeTable, totalnodecount);
            //printMatrix(SpTreeAdjMatrix, "first spanning tree");
            adjacencyMatrixOfNetwork = createAdjMatrix(currentnetwork, nodeList, edgeTable, totalnodecount, edgeWeightAttribute);
            //double[][] adjMatrixRemSPTree = new double[totalnodecount][totalnodecount];
            for (int i = 0; i < totalnodecount; i++) {
                for (int j = 0; j < totalnodecount; j++) {
                    if (SpTreeAdjMatrix[i][j] > 0) {
                        adjacencyMatrixOfNetwork[i][j] = 0.0;
                        adjacencyMatrixOfNetwork[j][i] = 0.0;
                        //adjMatrixRemSPTree[i][j] = adjacencyMatrixOfNetwork[i][j] - SpTreeAdjMatrix[i][j];
                    }

                }
            }
            //Matrix(adjacencyMatrixOfNetwork, "sptree removed matrix");
            SpTreeAdjMatrix = createSpTreeAdjMatrix(adjacencyMatrixOfNetwork, nodeList, edgeTable, totalnodecount);
            createNetwork(SpTreeAdjMatrix, nodeList, nodeTable, totalnodecount);
        }

    }

    public static double[][] createAdjMatrix(CyNetwork currentnetwork, List<CyNode> nodeList, CyTable edgeTable, int totalnodecount, String edgeWeightAttribute) {
        //make an adjacencymatrix for the current network
        double[][] adjacencyMatrixOfNetwork = new double[totalnodecount][totalnodecount];
        int k = 0;
        for (CyNode root : nodeList) {
            List<CyNode> neighbors = currentnetwork.getNeighborList(root, CyEdge.Type.OUTGOING);
            for (CyNode neighbor : neighbors) {
                List<CyEdge> edges = currentnetwork.getConnectingEdgeList(root, neighbor, CyEdge.Type.DIRECTED);
                if (edges.size() > 0) {
                    CyRow row = edgeTable.getRow(edges.get(0).getSUID());
                    try {
                        adjacencyMatrixOfNetwork[k][nodeList.indexOf(neighbor)] = Double.parseDouble(row.get(edgeWeightAttribute, String.class));
                    } catch (NumberFormatException ex) {
                    }
                }
            }
            k++;
        }
        return adjacencyMatrixOfNetwork;
    }

    public double[][] createSpTreeAdjMatrix(double[][] adjacencyMatrixOfNetwork, List<CyNode> nodeList, CyTable edgeTable, int totalnodecount) {
        // run kruskals algo
        double[][] adjacencyMatrixOfNetworkRelpica = new double[totalnodecount][totalnodecount];
        for (int i = 0; i < totalnodecount; i++) {
            for (int j = 0; j < totalnodecount; j++) {
                adjacencyMatrixOfNetworkRelpica[i][j] = adjacencyMatrixOfNetwork[i][j];
            }
        }
        Graph myGraph = new Graph(totalnodecount);
        int noOfEdges = 0;
        if (isMinimum) {
            while (noOfEdges < totalnodecount - 1) {
                Graph myGraphReplica = new Graph(myGraph);

                double minWeight = 32.0;
                int vertexi = 0;
                int vertexj = 0;
                for (int i = 0; i < totalnodecount; i++) {
                    for (int j = 0; j < totalnodecount; j++) {
                        if (adjacencyMatrixOfNetwork[i][j] != 0) {
                            if (adjacencyMatrixOfNetwork[i][j] < minWeight) {
                                minWeight = adjacencyMatrixOfNetwork[i][j];
                                vertexi = i;
                                vertexj = j;
                            }
                        }
                    }
                }
                myGraph.addEdge(vertexi, vertexj);
                adjacencyMatrixOfNetwork[vertexi][vertexj] = 0.0;
                if (minWeight == 32.0 && vertexi == 0 && vertexj == 0) {
                    break;
                }
                Cycle myGraphCycle = new Cycle(myGraph);
                if (myGraphCycle.hasCycle()) {
                    myGraph = myGraphReplica;
                } else {
                    noOfEdges++;
                }
            }
        } else {
            while (noOfEdges < totalnodecount - 1) {
                Graph myGraphReplica = new Graph(myGraph);

                double maxWeight = 0;
                int vertexi = 0;
                int vertexj = 0;
                for (int i = 0; i < totalnodecount; i++) {
                    for (int j = 0; j < totalnodecount; j++) {
                        if (adjacencyMatrixOfNetwork[i][j] > maxWeight) {
                            maxWeight = adjacencyMatrixOfNetwork[i][j];
                            vertexi = i;
                            vertexj = j;
                        }
                    }
                }
                myGraph.addEdge(vertexi, vertexj);
                adjacencyMatrixOfNetwork[vertexi][vertexj] = 0.0;
                //System.out.println("Max weight is "+maxWeight+" at i="+vertexi+", j="+vertexj);
                if (maxWeight == 0 && vertexi == 0 && vertexj == 0) {
                    break;
                }
                Cycle myGraphCycle = new Cycle(myGraph);
                if (myGraphCycle.hasCycle()) {
                    myGraph = myGraphReplica;
                } else {
                    noOfEdges++;
                }
            }
        }
        //make a new adjacency matrix for the output tree
        adjacencyMatrixOfNetwork = adjacencyMatrixOfNetworkRelpica;
        double[][] adjacencyMatrixOfNewNetwork = new double[totalnodecount][totalnodecount];
        for (int i = 0; i < totalnodecount; i++) {
            for (Integer j : myGraph.adj(i)) {
                adjacencyMatrixOfNewNetwork[i][j] = adjacencyMatrixOfNetwork[i][j];
            }
        }
        for (int i = 0; i < totalnodecount; i++) {
            for (int j = i + 1; j < totalnodecount; j++) {
                if (adjacencyMatrixOfNewNetwork[i][j] > adjacencyMatrixOfNewNetwork[j][i]) {
                    adjacencyMatrixOfNewNetwork[j][i] = 0.0;
                } else {
                    adjacencyMatrixOfNewNetwork[i][j] = 0.0;
                }
            }
        }
        //printMatrix(adjacencyMatrixOfNewNetwork, "inside");
        return adjacencyMatrixOfNewNetwork;
    }

    public void createNetwork(double[][] adjacencyMatrixOfNewNetwork, List<CyNode> nodeList, CyTable nodeTable, int totalnodecount) {
        CyNetwork SpanningTree;
        // To get a reference of CyNetworkFactory at CyActivator class of the App
        CyNetworkFactory networkFactory = CyActivator.networkFactory;
        // Create a new network
        SpanningTree = networkFactory.createNetwork();

        // Set name for network
        SpanningTree.getRow(SpanningTree).set(CyNetwork.NAME, "Spanning Tree");

        // Add nodes to the network
        List<CyNode> nodesInNewNetwork = new ArrayList<CyNode>(totalnodecount);
        for (int i = 0; i < nodeList.size(); i++) {
            nodesInNewNetwork.add(SpanningTree.addNode());
        }
        // Set name for new nodes
        for (int i = 0; i < nodeList.size(); i++) {
            SpanningTree.getRow(nodesInNewNetwork.get(i)).set(CyNetwork.NAME, nodeTable.getRow(nodeList.get(i).getSUID()).get(CyNetwork.NAME, String.class));
        }
        //add edges
        for (int i = 0; i < totalnodecount; i++) {
            for (int j = i + 1; j < totalnodecount; j++) {
                double maxi = Math.max(adjacencyMatrixOfNewNetwork[i][j], adjacencyMatrixOfNewNetwork[j][i]);
                if (maxi > 0.0) {
                    CyEdge root = SpanningTree.addEdge(nodesInNewNetwork.get(i), nodesInNewNetwork.get(j), true);
                    CyRow row = SpanningTree.getDefaultEdgeTable().getRow(root.getSUID());
                    row.set(edgeWeightAttribute, "" + maxi);
                }
            }
        }

        // Add the network to Cytoscape
        CyNetworkManager networkManager = CyActivator.networkManager;
        networkManager.addNetwork(SpanningTree);

        //Add view to cyto
        //CyNetworkView myView = CyActivator.networkViewFactory.createNetworkView(SpanningTree);
        //CyActivator.networkViewManager.addNetworkView(myView);
    }

    public static void printMatrix(double[][] matrix, String name) {
        System.out.println("The matrix " + name + " is :");
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    public static void printMatrixInt(int[][] matrix, String name) {
        System.out.println("The matrix " + name + " is :");
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }
    public void printGraph(Graph g, int n, String name) {
        System.out.println("The Graph " + name + " is :");
        for (int i = 0; i < n; i++) {
            System.out.print(i);
            for (Integer j : g.adj(i)) {
                System.out.print("--" + j);
            }
            System.out.println();
        }
    }
}
