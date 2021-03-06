/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cytoscape.spanningtree.internal.apsp;

import java.util.ArrayList;
import java.util.List;
import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.model.*;
import org.cytoscape.spanningtree.internal.CyActivator;
import org.cytoscape.spanningtree.internal.SpanningTreeStartMenu;
import org.cytoscape.view.model.CyNetworkView;

/**
 *
 * @author smd.faizan@gmail.com
 */
public class APSPalgoThread extends Thread {
    public boolean stop;
    CyApplicationManager cyApplicationManager;
    CyNetworkView currentnetworkview;
    CyNetwork currentnetwork;
    SpanningTreeStartMenu menu;

    public APSPalgoThread(CyApplicationManager cyApp, SpanningTreeStartMenu menu) {
        this.cyApplicationManager = cyApp;
        this.menu = menu;
    }

    @Override
    public void run() {
        stop = false;
        long lStartTime = System.currentTimeMillis();
        System.out.println("Start time for APSP algo: " + lStartTime + " milli seconds");
    
        currentnetworkview = cyApplicationManager.getCurrentNetworkView();
        currentnetwork = currentnetworkview.getModel();
        List<CyNode> nodeList = currentnetwork.getNodeList();
        int totalnodecount = nodeList.size();
        CyTable nodeTable = currentnetwork.getDefaultNodeTable();
        CyTable edgeTable = currentnetwork.getDefaultEdgeTable();

        double[][] adjacencyMatrix = createAdjMatrix(currentnetwork, nodeList, edgeTable, totalnodecount);
        if(adjacencyMatrix == null)
            return ;
        //SpanningTreeThread.printMatrix(adjacencyMatrix, "adjacency matrix");
        double[][] D = adjacencyMatrix.clone();
        int[][] Pi = runFloydAlgo(D);
        if(Pi == null)
            return;
        createNetwork(D, Pi, adjacencyMatrix, nodeList, nodeTable);

	long lEndTime = System.currentTimeMillis();
	long difference = lEndTime - lStartTime;
        System.out.println("End time for APSP algo: " + lEndTime + " milli seconds");
	System.out.println("Execution time for APSP algo: " + difference + " milli seconds");
        menu.setStatusLabel("All pair Shortest Paths graph created in Network panel");
    }

    public int[][] runFloydAlgo(double[][] D) {
        int totalnodecount = D.length;

        int[][] Pi = new int[totalnodecount][totalnodecount];
        for (int i = 0; i < totalnodecount; i++) {
            for (int j = 0; j < totalnodecount; j++) {
                if(stop)
                    return null;
                if (i == j || D[i][j] == 31999.0) {
                    Pi[i][j] = -1;
                } else {
                    Pi[i][j] = i;
                }
            }
        }
        for (int k = 0; k < totalnodecount; k++) {
            for (int i = 0; i < totalnodecount; i++) {
                for (int j = 0; j < totalnodecount; j++) {
                    if(stop)
                        return null;
                    if (i != j) {
                        if (D[i][j] > D[i][k] + D[k][j]) {
                            D[i][j] = D[i][k] + D[k][j];
                            Pi[i][j] = Pi[k][j];
                        }
                    }
                }
            }
        }
        return Pi;
    }

    public void createNetwork(double[][] D, int[][] Pi, double[][] adj, List<CyNode> nodeList, CyTable nodeTable) {
        int totalnodecount = D.length;
        CyNetwork Floyd;
        // To get a reference of CyNetworkFactory at CyActivator class of the App
        CyNetworkFactory networkFactory = CyActivator.networkFactory;
        // Create a new network
        Floyd = networkFactory.createNetwork();

        // Set name for network
        Floyd.getRow(Floyd).set(CyNetwork.NAME, "Floyd Warshall");

        // Add nodes to the network
        List<CyNode> nodesInNewNetwork = new ArrayList<CyNode>(totalnodecount);
        for (int i = 0; i < nodeList.size(); i++) {
            nodesInNewNetwork.add(Floyd.addNode());
        }
        // Set name for new nodes
        for (int i = 0; i < nodeList.size(); i++) {
            Floyd.getRow(nodesInNewNetwork.get(i)).set(CyNetwork.NAME, nodeTable.getRow(nodeList.get(i).getSUID()).get(CyNetwork.NAME, String.class));
        }
        //add edges
        for (int i = 0; i < totalnodecount; i++) {
            for (int j = i + 1; j < totalnodecount; j++) {
                if (D[i][j] != 31999.0) {
                    int k = Pi[i][j];
                    int to = j;
                    do {
                        // edge between k-->to
                        //System.out.println("k = "+k +", to = "+to);
                        List<CyEdge> edges = Floyd.getConnectingEdgeList(nodesInNewNetwork.get(k), nodesInNewNetwork.get(to), CyEdge.Type.ANY);
                        if (edges.isEmpty()) {
                            CyEdge root = Floyd.addEdge(nodesInNewNetwork.get(k), nodesInNewNetwork.get(to), true);
                            CyRow row = Floyd.getDefaultEdgeTable().getRow(root.getSUID());
                            row.set(CyEdge.INTERACTION, "" + adj[k][to]);
                        }
                        //edge added
                        to = k;
                        k = Pi[i][k];
                    } while (k != -1);
                }
            }
        }
        // Add the network to Cytoscape
        CyNetworkManager networkManager = CyActivator.networkManager;
        networkManager.addNetwork(Floyd);
        
        //Add view to cyto
//        CyNetworkView myView = CyActivator.networkViewFactory.createNetworkView(Floyd);
//        CyActivator.networkViewManager.addNetworkView(myView);
    }

    public double[][] createAdjMatrix(CyNetwork currentnetwork, List<CyNode> nodeList, CyTable edgeTable, int totalnodecount) {
        //make an adjacencymatrix for the current network
        double[][] adjacencyMatrixOfNetwork = new double[totalnodecount][totalnodecount];
        for (int i = 0; i < totalnodecount; i++) {
            for (int j = 0; j < totalnodecount; j++) {
                if(stop)
                    return null;
                adjacencyMatrixOfNetwork[i][j] = 31999.0;
            }
        }
        int k = 0;
        for (CyNode root : nodeList) {
            List<CyNode> neighbors = currentnetwork.getNeighborList(root, CyEdge.Type.ANY);
            for (CyNode neighbor : neighbors) {
                if(stop)
                    return null;
                List<CyEdge> edges = currentnetwork.getConnectingEdgeList(root, neighbor, CyEdge.Type.DIRECTED);
                if (edges.size() > 0) {
                    CyRow row = edgeTable.getRow(edges.get(0).getSUID());
                    try {
                        adjacencyMatrixOfNetwork[k][nodeList.indexOf(neighbor)] = Double.parseDouble(row.get(CyEdge.INTERACTION, String.class));
                    } catch (NumberFormatException ex) {
                    }
                }
            }
            k++;
        }
        return adjacencyMatrixOfNetwork;
    }
    
    public void end(){
        stop = true;
    }
    
}