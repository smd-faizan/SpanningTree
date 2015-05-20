package org.cytoscape.spanningtree.internal;
/*
 * SpanningTreeStartMenu.java
 *
 * Created on 22 Dec 2013, 18.38
 */

/**
 *
 * @author smd.faizan@gmail.com
 */
import org.cytoscape.spanningtree.internal.apsp.APSPalgoThread;
import java.awt.Component;
import javax.swing.Icon;
import javax.swing.JOptionPane;
import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.application.swing.CySwingApplication;
import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.application.swing.CytoPanelName;
import org.cytoscape.model.*;
import org.cytoscape.spanningtree.internal.spanningtree.SpanningTreeThread;
import org.cytoscape.view.model.CyNetworkView;

public class SpanningTreeStartMenu extends javax.swing.JPanel implements CytoPanelComponent {

    private SpanningTreeCore spanningtreecore;
    CyApplicationManager cyApplicationManager;
    CySwingApplication cyDesktopService;
    CyNetwork currentnetwork;
    CyNetworkView currentnetworkview;
    public CyActivator cyactivator;
    static String edgeWeightAttribute;
    public SpanningTreeThread spannigTreeThread;
    public APSPalgoThread algo;

    public SpanningTreeStartMenu(CyActivator cyactivator, SpanningTreeCore spanningtreecore) {
        initComponents();
        this.cyactivator = cyactivator;
        this.spanningtreecore = spanningtreecore;
        cyApplicationManager = spanningtreecore.getCyApplicationManager();
        cyDesktopService = spanningtreecore.getCyDesktopService();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        buttonGroup2 = new javax.swing.ButtonGroup();
        jScrollPane1 = new javax.swing.JScrollPane();
        jPanel5 = new javax.swing.JPanel();
        jPanel2 = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        spanningTreeButton = new javax.swing.JButton();
        alternateSpanningTreeButton = new javax.swing.JButton();
        maxRadioButton = new javax.swing.JRadioButton();
        minRadioButton = new javax.swing.JRadioButton();
        exitButton = new javax.swing.JButton();
        statusLabel = new javax.swing.JLabel();
        helpButton = new javax.swing.JButton();
        APSPalgoButton = new javax.swing.JButton();
        statusBar = new javax.swing.JProgressBar();
        stopButton = new javax.swing.JButton();
        jLabel2 = new javax.swing.JLabel();

        setBorder(javax.swing.BorderFactory.createEtchedBorder());
        setRequestFocusEnabled(false);

        jScrollPane1.setToolTipText("");
        jScrollPane1.setMaximumSize(new java.awt.Dimension(400, 32767));
        jScrollPane1.setPreferredSize(new java.awt.Dimension(300, 860));

        jPanel5.setBorder(new javax.swing.border.MatteBorder(null));

        jPanel2.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        jLabel1.setText("Click the buttons to make a spanning tree network");

        spanningTreeButton.setText("Create a Spanning Tree");
        spanningTreeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                spanningTreeButtonActionPerformed(evt);
            }
        });

        alternateSpanningTreeButton.setText("Create Second Spanning Tree");
        alternateSpanningTreeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                alternateSpanningTreeButtonActionPerformed(evt);
            }
        });

        buttonGroup2.add(maxRadioButton);
        maxRadioButton.setText("Maximal");

        buttonGroup2.add(minRadioButton);
        minRadioButton.setSelected(true);
        minRadioButton.setText("Minimal");

        exitButton.setText("Exit");
        exitButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exitButtonActionPerformed(evt);
            }
        });

        statusLabel.setText("Tree Making Status");

        helpButton.setText("Help");
        helpButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                helpButtonActionPerformed(evt);
            }
        });

        APSPalgoButton.setText("Create All-Pair-Shortest-Path graph");
        APSPalgoButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                APSPalgoButtonActionPerformed(evt);
            }
        });

        stopButton.setText("Stop");
        stopButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                stopButtonActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout jPanel2Layout = new org.jdesktop.layout.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .add(jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(jLabel1, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.TRAILING, jPanel2Layout.createSequentialGroup()
                        .add(jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, statusBar, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .add(APSPalgoButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .add(jPanel2Layout.createSequentialGroup()
                                .add(spanningTreeButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(alternateSpanningTreeButton)))
                        .add(4, 4, 4))
                    .add(jPanel2Layout.createSequentialGroup()
                        .add(jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(maxRadioButton)
                            .add(minRadioButton))
                        .add(0, 0, Short.MAX_VALUE))
                    .add(jPanel2Layout.createSequentialGroup()
                        .add(jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(statusLabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .add(jPanel2Layout.createSequentialGroup()
                                .add(helpButton, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 77, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                        .add(jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
                            .add(stopButton, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 82, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, exitButton, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 82, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                        .add(2, 2, 2)))
                .addContainerGap())
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(org.jdesktop.layout.GroupLayout.TRAILING, jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .add(jLabel1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 14, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                .add(minRadioButton, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 23, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(maxRadioButton)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                .add(jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(spanningTreeButton)
                    .add(alternateSpanningTreeButton))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                .add(APSPalgoButton)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, 18, Short.MAX_VALUE)
                .add(statusBar, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 23, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(statusLabel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 21, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(stopButton))
                .add(18, 18, 18)
                .add(jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(helpButton)
                    .add(exitButton))
                .addContainerGap())
        );

        jLabel2.setFont(new java.awt.Font("Dialog", 1, 16)); // NOI18N
        jLabel2.setForeground(new java.awt.Color(255, 0, 51));
        jLabel2.setText("Spanning tree Menu");

        org.jdesktop.layout.GroupLayout jPanel5Layout = new org.jdesktop.layout.GroupLayout(jPanel5);
        jPanel5.setLayout(jPanel5Layout);
        jPanel5Layout.setHorizontalGroup(
            jPanel5Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jPanel5Layout.createSequentialGroup()
                .add(22, 22, 22)
                .add(jPanel5Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(jLabel2)
                    .add(jPanel2, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .add(0, 812, Short.MAX_VALUE))
        );
        jPanel5Layout.setVerticalGroup(
            jPanel5Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jPanel5Layout.createSequentialGroup()
                .addContainerGap()
                .add(jLabel2)
                .add(18, 18, 18)
                .add(jPanel2, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(899, Short.MAX_VALUE))
        );

        jScrollPane1.setViewportView(jPanel5);

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(jScrollPane1, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 407, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(jScrollPane1, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1010, Short.MAX_VALUE)
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents

    private void spanningTreeButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_spanningTreeButtonActionPerformed
        currentnetworkview = cyApplicationManager.getCurrentNetworkView();
        currentnetwork = currentnetworkview.getModel();
        edgeWeightAttribute = inputEdgeAttributeAndValidate(currentnetwork.getDefaultEdgeTable());
        if(edgeWeightAttribute == null)
            return;
        statusLabel.setText("Making Spanning Tree...");
        spannigTreeThread = new SpanningTreeThread(currentnetwork, currentnetworkview, false, minRadioButton.isSelected(), edgeWeightAttribute, this);
        spannigTreeThread.start();
    }//GEN-LAST:event_spanningTreeButtonActionPerformed

    private void alternateSpanningTreeButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_alternateSpanningTreeButtonActionPerformed
        currentnetworkview = cyApplicationManager.getCurrentNetworkView();
        currentnetwork = currentnetworkview.getModel();
        edgeWeightAttribute = inputEdgeAttributeAndValidate(currentnetwork.getDefaultEdgeTable());
        if(edgeWeightAttribute == null)
            return;
        statusLabel.setText("Making Spanning Tree...");
        spannigTreeThread = new SpanningTreeThread(currentnetwork, currentnetworkview, true, minRadioButton.isSelected(), edgeWeightAttribute, this );
        spannigTreeThread.start();
    }//GEN-LAST:event_alternateSpanningTreeButtonActionPerformed

    private void exitButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exitButtonActionPerformed
        // TODO add your handling code here:
        spanningtreecore.closecore();
        spanningtreecore.closeSpanningTreeStartMenu();
    }//GEN-LAST:event_exitButtonActionPerformed

    private void helpButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_helpButtonActionPerformed
        // TODO add your handling code here:
        SpanningTreeHelp help = new SpanningTreeHelp();
        help.setText(1);
        help.setVisible(true);
    }//GEN-LAST:event_helpButtonActionPerformed

    private void APSPalgoButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_APSPalgoButtonActionPerformed
        currentnetworkview = cyApplicationManager.getCurrentNetworkView();
        currentnetwork = currentnetworkview.getModel();
        edgeWeightAttribute = inputEdgeAttributeAndValidate(currentnetwork.getDefaultEdgeTable());
        if(edgeWeightAttribute == null)
            return;
        statusLabel.setText("Running All pair Shortest path Algorithm...");
        algo = new APSPalgoThread(cyApplicationManager, this);
        algo.start();
    }//GEN-LAST:event_APSPalgoButtonActionPerformed

    private void stopButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_stopButtonActionPerformed
        // TODO add your handling code here:
        if(spannigTreeThread.isAlive()){
            spannigTreeThread.end();
            stopcalculus(null);
        }
        if(algo.isAlive()){
            algo.end();
            stopcalculus(null);
        }
    }//GEN-LAST:event_stopButtonActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton APSPalgoButton;
    private javax.swing.JButton alternateSpanningTreeButton;
    private javax.swing.ButtonGroup buttonGroup2;
    private javax.swing.JButton exitButton;
    private javax.swing.JButton helpButton;
    private static javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel5;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JRadioButton maxRadioButton;
    private javax.swing.JRadioButton minRadioButton;
    private javax.swing.JButton spanningTreeButton;
    private javax.swing.JProgressBar statusBar;
    private javax.swing.JLabel statusLabel;
    private javax.swing.JButton stopButton;
    // End of variables declaration//GEN-END:variables

    @Override
    public Component getComponent() {
        return this;
    }

    @Override
    public CytoPanelName getCytoPanelName() {
        return CytoPanelName.WEST;
    }

    @Override
    public String getTitle() {
        return "Spanning Tree";
    }

    @Override
    public Icon getIcon() {
        return null;
    }
    
    public void setStatusLabel(String status){
        this.statusLabel.setText(status);
    }
    
    public String inputEdgeAttributeAndValidate(CyTable edgeTable){
        //System.out.println("Waiting for user input of edge name..");
        edgeWeightAttribute = JOptionPane.showInputDialog(null, "Enter the Edge Attribute to be used as edge weight for the network (case sensitive).");
        if(edgeWeightAttribute == null || edgeWeightAttribute.equals("")){
            System.out.println("Trying to use 'weight' or distance' or 'CyEdge.INTERACTION' as edge attribute ");
            if(edgeTable.getColumn("weight") != null){
                System.out.println("using 'weight' as edge attribute.");
                return "weight";
            } else if(edgeTable.getColumn("distance") != null){
                System.out.println("using 'distance' as edge attribute.");
                return "distance";
            } else if(edgeTable.getColumn(CyEdge.INTERACTION) != null){
                System.out.println("using 'CyEdge.INTERACTION' as edge attribute.");
                return CyEdge.INTERACTION;
            } else {
                JOptionPane.showMessageDialog(null, " no column with name 'weight' or distance' or 'CyEdge.INTERACTION' exists. Exiting!", "Spanning Tree", JOptionPane.WARNING_MESSAGE);
                return null;
            }
        } else {
            if(edgeTable.getColumn(edgeWeightAttribute) != null){
                System.out.println("using "+edgeWeightAttribute+" as edge attribute.");
                return edgeWeightAttribute;
            } else{
                JOptionPane.showMessageDialog(null, " no column with name "+edgeWeightAttribute+" exists. Exiting!", "Spanning Tree", JOptionPane.WARNING_MESSAGE);
                return null;         
            }
        }
    }
    
    
    public void endOfComputation(String message) {
    	statusBar.setIndeterminate(false);
        if(message == null)
            statusLabel.setText("Calculation Done.");
        else
            statusLabel.setText(message);
        spanningTreeButton.setEnabled(true);
    }

    public void stopcalculus(String message) {
        statusBar.setIndeterminate(false);
        if(message == null)
            statusLabel.setText("Interrupted by the user, click start to repeat");
        else
            statusLabel.setText(message);
        spanningTreeButton.setEnabled(true);
    }

    public void calculatingresult(String message) {
        statusBar.setIndeterminate(true);
        statusBar.setVisible(true);
        if(message == null)
            statusLabel.setText("Working ...");
        else
            statusLabel.setText(message);
        spanningTreeButton.setEnabled(false);
    }
}
