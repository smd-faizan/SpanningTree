/*
 * SpanningTreeHelp.java
 *
 * Created on 22 Dec 2013, 21.19
 */
package org.cytoscape.spanningtree.internal;

/**
 *
 * @author smd.faizan@gmail.com
 */
public class SpanningTreeHelp extends javax.swing.JFrame {

    private String helpString;

    /**
     * Creates new form SpanningTreeHelp
     */
    public SpanningTreeHelp() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {
        jScrollPane1 = new javax.swing.JScrollPane();
        jTextArea1 = new javax.swing.JTextArea();
        setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        setForeground(new java.awt.Color(255, 255, 51));
        jTextArea1.setColumns(20);
        jTextArea1.setRows(5);
        jScrollPane1.setViewportView(jTextArea1);
        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
                layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(layout.createSequentialGroup().addContainerGap().add(jScrollPane1, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 518, Short.MAX_VALUE).addContainerGap()));
        layout.setVerticalGroup(
                layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(layout.createSequentialGroup().addContainerGap().add(jScrollPane1, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 520, Short.MAX_VALUE).addContainerGap()));
        pack();
    }// </editor-fold>//GEN-END:initComponents

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        java.awt.EventQueue.invokeLater(new Runnable() {

            public void run() {
                new SpanningTreeHelp().setVisible(true);
            }
        });
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JTextArea jTextArea1;
    // End of variables declaration//GEN-END:variables

    public void setText(int buttonNumber) {
        String secondSpanningTree = 
                "Second Spanning Tree:\n\n\n"
                + "Spanning tree which does not include any of the edges form the first \n"
                + "spanning tree. It is calculated as below.\n"
                + "First create a graph by removing all the edges of first spanning tree \n"
                + "from the original graph. Now calculate the spanning tree for this new \n"
                + "graph.\n\n"
                ;
        String APSP = 
                "All-Pair-Shortest-Path  graph:\n\n\n"
                + "All-Pair-Shortest-Path also known as Floyd Warshall Algorithm. Deatils \n"
                + "can be found in the wiki http://en.wikipedia.org/wiki/Floyd-Warshall_algorithm\n\n"
                ;
        helpString =
                "Tree: \n\n"
                + "Tree is a minimally connected graph. Specifically in graph theory,\n"
                + "a tree is an undirected graph in which any two vertices are connected\n"
                + "by exactly one simple path. In other words, any connected graph without\n"
                + "simple cycles is a tree.\n\n"
                + "Spanning Tree :\n\n\n"
                + "Given a graphA which is not a tree, we can make graphA a tree by eliminating\n"
                + "some edges and thus we can be able to construct a tree out of the graphA.\n"
                + "This problem becomes interesting when edges have weights. Now, removing\n"
                + "different edges in graphA results in different trees. In a minimum spanning tree,\n"
                + "one would ideally remove edges with high weights (or) low weights depending upon\n"
                + "the problem.When it comes to simulation of transportation network where\n"
                + "edges represent the distance between cities, one would like to remove edges\n"
                + "with high weights(distances) to form spanning tree. When it comes to simulation\n"
                + "of wireless sensor network where edges represent the edge connectivity, one\n"
                + "would like to remove edges which have less edge connectivity to ensure reliability\n"
                + "constraints.\n\n"
                + secondSpanningTree
                + APSP;
        this.setTitle("Spanning Tree Help : ");
        
        jTextArea1.setText(helpString);
        jTextArea1.setCaretPosition(0);
    }
}
