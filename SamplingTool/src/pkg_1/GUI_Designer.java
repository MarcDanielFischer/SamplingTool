/*
 *      
 *      
 */

package pkg_1;

import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.table.DefaultTableModel;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.swing.JScrollPane;
import javax.swing.JList;
import javax.swing.JTable;
import javax.swing.JTextField;


@SuppressWarnings("serial")
public class GUI_Designer extends JFrame {
	
	// GUI components
	private JPanel contentPane;
	private JFileChooser openFileDialog;
	private JButton btnBrowse;
	private JButton btnRunSampling; 
	private JButton btnClose;
	private JTextField textField_Start_X;
	private JTextField textField_Start_Y;
	private JTextField textField_NumSubPlotsPerCluster;
	private JTextField textFieldDistBetweenSubPlots;
	private JTextField textField_H_numPlotsVertical;
	private JTextField textField_DistX;
	private JTextField textField_DistY;
	private JLabel lblSelectedFile;
	private JLabel lblSubplotsVerticalLine;
	private JLabel lblClusterDesign;
	private JLabel lblDistanceBetweenSubplots;
	private JLabel lblSubplotsPerCluster; 
	private JLabel lblChooseStrataFor;
	private JLabel lblSamplingDesign;
	private JLabel lblStart_Y; 
	private JLabel lblNewLabel; 
	private JLabel lbl_SelectColumn;
	private JLabel lblDist_x;
	private JLabel lblStartingPoint;
	private JLabel lblNewLabel_3;
	private JLabel lblStart_X;
	private JLabel lblDist_y;
	private JComboBox comboBox_StartingPoint;
	private JComboBox comboBox_ClusterSampling;
	private JComboBox comboBox_Columns;
	private JComboBox comboBox_ClusterDesign;
	private JComboBox comboBox_SamplingDesign;
	private JScrollPane scrollPaneStrataList;
	private JList strataNamesList;
	private DefaultTableModel model;
	private JTable table;
	private String[] strataNamesArray;
	
	// Sampling parameters
	// alle nicht in jedem Fall ben�tigten Parameter habe ich schon hier initialisiert,
	// weil (momentan) immer alle Params an die Methode runSampling() �bergeben werden
	// und sie daf�r initialisiert sein m�ssen
	File inputFile;
	String sampleColumn;
	int numStrata;
	String[] selectedStrata;
	int samplingDesign;
	static final int SIMPLE_RANDOM_SAMPLING = 1; // constants to specify samplingDesign
	int[] numPlotsToBeSampled = null; // only for option SIMPLE_RANDOM_SAMPLING 
	static final int SYSTEMATIC_SAMPLING = 2;
	int gridDistX = 0; // only for option SYSTEMATIC_SAMPLING 
	int gridDistY = 0;// only for option SYSTEMATIC_SAMPLING 
	int startingPoint = 0; // only for option SYSTEMATIC_SAMPLING 
	static final int STARTING_POINT_RANDOM = 1; // constants to specify startingPoint
	static final int STARTING_POINT_SPECIFIED = 2;
	int startX = 0; // only for option STARTING_POINT_SPECIFIED 
	int startY = 0; // only for option STARTING_POINT_SPECIFIED 
	int clusterSampling;
	static final int CLUSTER_SAMPLING_NO = 1; // constants to specify clusterSampling
	static final int CLUSTER_SAMPLING_YES = 2;
	int clusterShape = 0; // only for clusterSampling == YES
	static final int I_SHAPE = 1; // constants to specify clusterShape
	static final int L_SHAPE = 2;
	static final int H_SHAPE = 3;
	static final int SQUARE_SHAPE = 4;
	int numSubPlotsinHVerticalLine = 0; // only for clusterShape == H_SHAPE
	int numClusterSubPlots = 0; // only for clusterSampling == YES
	int distBetweenSubPlots = 0; // only for clusterSampling == YES
	private JTextField textField_H_numPlotsHorizontal;
	
	
	
	
	
	
	

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		// dieses ganze Ged�ns mit EventQueue, invokeLater() und run() ist f�r die Programmausf�hrung durch Threads gedacht
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					GUI_Designer frame = new GUI_Designer();
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 */
	public GUI_Designer() {
		// Application Window Settings
		setTitle("Arbonaut Spatial Sampling Tool");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 644, 902);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		// Cancel Button
		btnClose = new JButton("Close");
		btnClose.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				System.exit(0); // evtl. checken, ob das die Option der Wahl ist
			}
		});
		btnClose.setBounds(496, 800, 122, 23);
		contentPane.add(btnClose);
		
		
		
		// "Browse" button
		btnBrowse = new JButton("Browse");
		btnBrowse.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				openFileDialog = new JFileChooser();
				if(openFileDialog.showOpenDialog(null) == JFileChooser.APPROVE_OPTION){
					inputFile = openFileDialog.getSelectedFile();
					// write input file path to Label (just for visualization purposes)
					lblSelectedFile.setText(inputFile.toString());
					
					// write input SHP column names to comboBox_StratField immediately after File is selected
					try{
						// clear comboBox_StratField entries first if selected file has been changed
						comboBox_Columns.removeAllItems();
						String[] colNames = SamplingFunctionalityMethods.getSHPColNames(inputFile);
						// add SHP column names to comboBox_StratField dropdown menu
						for(String colName : colNames){
							// TODO sort items alphabetically when adding them to comboBox_StratField
							comboBox_Columns.addItem(colName);
						}
					}catch(Exception e){
						JOptionPane.showMessageDialog(null, "Some error ocurred, possibly you have not selected a Shapefile" + e.toString());
					}
				}
			}
		});
		btnBrowse.setBounds(10, 25, 89, 23);
		contentPane.add(btnBrowse);
		
		
		
		// Combo Box "Columns"
		comboBox_Columns = new JComboBox();
		comboBox_Columns.setBounds(10, 145, 200, 30);
		contentPane.add(comboBox_Columns);
		// when a column is selected, all the values in this column will be read automatically and filled as 
        // items into the "Strata" JList
		comboBox_Columns.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				// String actionCommand = e.getActionCommand(); // actionCommand "comboBoxChanged" kommt sowohl bei addItem() als auch bei Item-Auswahl
				JComboBox eventSource = (JComboBox)e.getSource();
				String selectedItem = (String) eventSource.getSelectedItem();
				// Workaround: beim Bef�llen der ComboBox wird sofort ein ActionEvent ausgel�st, obwohl das erst passieren soll, wenn der User ein Item selected
				//  --> also �ber if-Block filtern
				if(selectedItem != "the_geom" && selectedItem != null){ // die Werte f�r Spalte the_geom sollen nicht ausgelesen werden
					try{
						// read values from specified column and write them to JList
						ArrayList<String> strataNames = SamplingFunctionalityMethods.getStrataColumnValues(inputFile, selectedItem);
						// sort strataNames alphabetically 
						Collections.sort(strataNames);
						strataNamesList = new JList(strataNames.toArray());
						scrollPaneStrataList.setViewportView(strataNamesList);


					}catch(Exception ex){
						JOptionPane.showMessageDialog(null, ex.toString());
					}
				}
			}
		});
		
        // Combo Box "Sampling Design"
		comboBox_SamplingDesign = new JComboBox();
		comboBox_SamplingDesign.setBounds(10, 401, 200, 30);
		comboBox_SamplingDesign.addItem(null); // damit die ComboBox am Anfang leer ist --> geht vl noch sch�ner
		comboBox_SamplingDesign.addItem("Simple Random Sampling");
		comboBox_SamplingDesign.addItem("Systematic Sampling");
		contentPane.add(comboBox_SamplingDesign);
		comboBox_SamplingDesign.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				/*
				 * Sampling-Logik in GUI umsetzen: bedingte Anzeige von Komponenten:
				 * 
				 * - TextField "Number of Plots" soll nur angezeigt werden bei Simple Random Sampling
				 * - Eingabeparameter "Dist_X", "Dist_Y" und "Starting Point" sollen nur bei systematic Sampling angezeigt werden
				 * 
				 * --> zum Hinzuf�gen muss ich repaint() aufrufen, geht aber nicht von hier aus --> Methode gekapselt und ausgelagert
				 * --> wenn ich die Sampling-Methode �ndere, bleiben die Kompnenten der anderen Option erhalten --> wie gehen die wieder weg?
				 */
				
				if((String)comboBox_SamplingDesign.getSelectedItem() == null){ // direkter Test auf "null" --> sieht irgendwie falsch aus
					disableSamplingDesignControls();
				}
				if((String)comboBox_SamplingDesign.getSelectedItem() == "Simple Random Sampling"){
					// disable previously enabled Sampling Design Controls, if systematic Sampling option had been selected before
					disableSamplingDesignControls();
				}
				
				if((String)comboBox_SamplingDesign.getSelectedItem() == "Systematic Sampling"){
					// disable previously enabled Sampling Design Controls, if Simple Random Sampling option had been selected before
					disableSamplingDesignControls();
					// selectively enable controls
					textField_DistX.setEnabled(true);
					lblDist_x.setEnabled(true);
					textField_DistY.setEnabled(true);
					lblDist_y.setEnabled(true);
					comboBox_StartingPoint.setEnabled(true);
					lblStartingPoint.setEnabled(true);
					
				
					}
					

				
			}
		});
		
		// Combo Box "Starting Point"
		comboBox_StartingPoint = new JComboBox();
		comboBox_StartingPoint.setEnabled(false);
		comboBox_StartingPoint.setBounds(10, 642, 183, 30);
		comboBox_StartingPoint.addItem(null);
		comboBox_StartingPoint.addItem("Random");
		comboBox_StartingPoint.addItem("Specified");
		contentPane.add(comboBox_StartingPoint);
		comboBox_StartingPoint.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				if((String) comboBox_StartingPoint.getSelectedItem() == "Specified"){
					// enable controls for Systematic Sampling Specified Starting Point X and Y coordinate Selection
					// --> die folgenden 4 Zeilen k�nnte ich evtl. in eine MEthode auslagern, weil sie doppelt vorkommen
					lblStart_X.setEnabled(true);
					lblStart_Y.setEnabled(true);
					textField_Start_X.setEnabled(true);
					textField_Start_Y.setEnabled(true);
				}
				// man kann direkt auf null testen: var == null --> stimmt das so? es funktioniert jedenfalls
				if((String) comboBox_StartingPoint.getSelectedItem() == "Random" || (String) comboBox_StartingPoint.getSelectedItem() == null){
					// disable controls for Systematic Sampling Specified Starting Point X and Y coordinate Selection
					lblStart_X.setEnabled(false);
					lblStart_Y.setEnabled(false);
					textField_Start_X.setEnabled(false);
					textField_Start_Y.setEnabled(false);
				}
				
			}
		});
		
		
		// ComboBox "Cluster Sampling"
		comboBox_ClusterSampling = new JComboBox();
		comboBox_ClusterSampling.setBounds(221, 401, 200, 30);
		comboBox_ClusterSampling.addItem(null);
		comboBox_ClusterSampling.addItem("Yes");
		comboBox_ClusterSampling.addItem("No");
		contentPane.add(comboBox_ClusterSampling);
		comboBox_ClusterSampling.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				// if ClusterSampling==Yes then enable controls
				if(comboBox_ClusterSampling.getSelectedItem() == "Yes"){
					// disable all previously enabled Cluster Controls to avoid confusion (accidentally enabled controls)
					disableClusterControls();
					// enable controls
					lblClusterDesign.setEnabled(true);
					comboBox_ClusterDesign.setEnabled(true);
					lblSubplotsPerCluster.setEnabled(true);
					textField_NumSubPlotsPerCluster.setEnabled(true);
					lblDistanceBetweenSubplots.setEnabled(true);
					textFieldDistBetweenSubPlots.setEnabled(true);
					
				}else{ // if ClusterSampling==No OR none selected then disable all controls
					disableClusterControls();
				}
			}
		});
		
		
		// ComboBox "Cluster Design"
		comboBox_ClusterDesign = new JComboBox();
		comboBox_ClusterDesign.setEnabled(false);
		comboBox_ClusterDesign.setBounds(221, 460, 200, 30);
		// add selectable cluster shape entries: I,L,H,Square
		comboBox_ClusterDesign.addItem(null);
		comboBox_ClusterDesign.addItem("I");
		comboBox_ClusterDesign.addItem("L");
		comboBox_ClusterDesign.addItem("H");
		comboBox_ClusterDesign.addItem("Square");
		contentPane.add(comboBox_ClusterDesign);
		comboBox_ClusterDesign.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				// if Cluster Shape == H then enable "number of Plots in one vertical line" Control else disable said control
				if(comboBox_ClusterDesign.getSelectedItem()=="H"){
					lblSubplotsVerticalLine.setEnabled(true);
					textField_H_numPlotsVertical.setEnabled(true);
				}else{
					lblSubplotsVerticalLine.setEnabled(false);
					textField_H_numPlotsVertical.setEnabled(false);
				}
			}
		});
		
		
		// "Run Sampling" Button 
		btnRunSampling = new JButton("Run Sampling");
		btnRunSampling.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				// alle eingetragenen Samplingparameter zusammentragen und auf die am Klassenanfang deklarierten Variablen schreiben
				// TODO besseres Verfahren �berlegen, die ganzen Samplingparameter zu �bergeben (zu viele f�r eine einzelne Methode) --> Hashmap, eigenes Objekt nur f�r die Parameter?
				
				//----------------------------------------
				// Start Samplingparameter
				// --> Erfassung der Parameter enth�lt Samplinglogik
				sampleColumn = (String)comboBox_Columns.getSelectedItem();
				numStrata = model.getRowCount();
				selectedStrata = new String[numStrata];
				// read stratum names from JTable
				for(int i = 0; i < numStrata; i++){ // iterate over rows
					selectedStrata[i] = (String) model.getValueAt(i, 0); // Spalte "Name"
				}

				// Sampling Design params
				if((String)comboBox_SamplingDesign.getSelectedItem()=="Simple Random Sampling"){
					samplingDesign = SIMPLE_RANDOM_SAMPLING;
					// Number of plots already specified in Array "numPlotsToBeSampled"
					numPlotsToBeSampled = new int[numStrata];
					for(int i = 0; i < numStrata; i++){ // iterate over rows
						numPlotsToBeSampled[i] = Integer.parseInt((String)model.getValueAt(i, 1)); // Spalte "numPlots"
					}
				}else if ((String)comboBox_SamplingDesign.getSelectedItem()=="Systematic Sampling"){
					samplingDesign = SYSTEMATIC_SAMPLING;
					gridDistX = Integer.parseInt(textField_DistX.getText()); 
					gridDistY = Integer.parseInt(textField_DistY.getText()); 
					// Starting Point 
					if((String)comboBox_StartingPoint.getSelectedItem() == "Specified"){
						startingPoint = STARTING_POINT_SPECIFIED;
						startX = Integer.parseInt(textField_Start_X.getText());
						startY = Integer.parseInt(textField_Start_Y.getText());
					} else{
						startingPoint = STARTING_POINT_RANDOM; // use STARTING_POINT_RANDOM as default unless otherwise specified
					}
				}
				
				
				// Cluster Sampling params
				if((String)comboBox_ClusterSampling.getSelectedItem()=="Yes"){
					clusterSampling = CLUSTER_SAMPLING_YES;
					// all following params are only needed if clusterSampling == YES,
					// hence they are treated inside this if statement
					// Cluster Shape
					switch((String)comboBox_ClusterDesign.getSelectedItem()){
					case "I": clusterShape = I_SHAPE;
					break;
					case "L": clusterShape = L_SHAPE;
					break;
					case "H": clusterShape = H_SHAPE;
					break;
					case "Square": clusterShape = SQUARE_SHAPE;
					break;
					}
					// Number of Plots in a vertical line of the H-Shape
					if(clusterShape == H_SHAPE){
						numSubPlotsinHVerticalLine = Integer.parseInt(textField_H_numPlotsVertical.getText());
					}
					// number of sub-plots per cluster
					numClusterSubPlots = Integer.parseInt(textField_NumSubPlotsPerCluster.getText());
					// Distance between Cluster SubPlots
					distBetweenSubPlots = Integer.parseInt(textFieldDistBetweenSubPlots.getText());

				}else clusterSampling = CLUSTER_SAMPLING_NO; // use CLUSTER_SAMPLING_NO as default option, unless CLUSTER_SAMPLING_YES is explicitly specified
				// Ende Samplingparameter
				//----------------------------------------

				// Methode erst hinterher aufrufen mit gesammelten Sachen drin --> Methodensignatur �ndern
				// TODO Methode vl mit weniger Parameter hinbekommen --> evtl. Params als eigene Objektklasse
				boolean samplingSuccessful = SamplingFunctionalityMethods.runSampling(inputFile, sampleColumn, numStrata, selectedStrata, samplingDesign, numPlotsToBeSampled, gridDistX, gridDistY, startingPoint, startX, startY, clusterSampling, clusterShape, numSubPlotsinHVerticalLine, numClusterSubPlots, distBetweenSubPlots ); 
				if(samplingSuccessful){
					JOptionPane.showMessageDialog(null, "output file successfully written");
				}else{
					JOptionPane.showMessageDialog(null, "some error occurred");
				}



				}
			});
		btnRunSampling.setBounds(4, 800, 122, 23);
		contentPane.add(btnRunSampling);

		// TextField "Grid Distance between points in X direction"
		textField_DistX = new JTextField();
		textField_DistX.setEnabled(false);
		textField_DistX.setBounds(10, 536, 86, 20);
		contentPane.add(textField_DistX);
		textField_DistX.setColumns(10);

		// TextField "Grid Distance between points in Y direction"
		textField_DistY = new JTextField();
		textField_DistY.setEnabled(false);
		textField_DistY.setBounds(10, 590, 86, 20);
		contentPane.add(textField_DistY);
		textField_DistY.setColumns(10);
				
		// TextField "Start X"
		textField_Start_X = new JTextField();
		textField_Start_X.setEnabled(false);
		textField_Start_X.setBounds(10, 703, 86, 20);
		contentPane.add(textField_Start_X);
		textField_Start_X.setColumns(10);

		// TextField "Start Y"
		textField_Start_Y = new JTextField();
		textField_Start_Y.setEnabled(false);
		textField_Start_Y.setColumns(10);
		textField_Start_Y.setBounds(10, 759, 86, 20);
		contentPane.add(textField_Start_Y);
		
		// TextField "NumSubPlotsPerCluster"
		textField_NumSubPlotsPerCluster = new JTextField();
		textField_NumSubPlotsPerCluster.setEnabled(false);
		textField_NumSubPlotsPerCluster.setBounds(221, 536, 86, 20);
		contentPane.add(textField_NumSubPlotsPerCluster);
		textField_NumSubPlotsPerCluster.setColumns(10);
		
		// TextField "DistBetweenSubPlots"
		textFieldDistBetweenSubPlots = new JTextField();
		textFieldDistBetweenSubPlots.setEnabled(false);
		textFieldDistBetweenSubPlots.setBounds(221, 590, 86, 20);
		contentPane.add(textFieldDistBetweenSubPlots);
		textFieldDistBetweenSubPlots.setColumns(10);
		
		// TextField "H_numPlotsVertical"
		textField_H_numPlotsVertical = new JTextField();
		textField_H_numPlotsVertical.setEnabled(false);
		textField_H_numPlotsVertical.setBounds(221, 647, 86, 20);
		contentPane.add(textField_H_numPlotsVertical);
		textField_H_numPlotsVertical.setColumns(10);
		// End TextFields
		//-----------------------------------------------------------------------------------------
		
		
		
		
		//-----------------------------------------------------------------------------------------
		// Labels
		lblSelectedFile = new JLabel("No File selected yet");
		lblSelectedFile.setBounds(92, 47, 448, 30);
		contentPane.add(lblSelectedFile);
		
		lblSubplotsVerticalLine = new JLabel("Sub-plots per vertical line (only H Clusters):");
		lblSubplotsVerticalLine.setEnabled(false);
		lblSubplotsVerticalLine.setBounds(221, 602, 300, 50);
		contentPane.add(lblSubplotsVerticalLine);
		
		lblClusterDesign = new JLabel("Cluster Design:");
		lblClusterDesign.setEnabled(false);
		lblClusterDesign.setBounds(217, 416, 200, 50);
		contentPane.add(lblClusterDesign);

		lblDistanceBetweenSubplots = new JLabel("Distance between sub-plots");
		lblDistanceBetweenSubplots.setEnabled(false);
		lblDistanceBetweenSubplots.setBounds(221, 554, 200, 50);
		contentPane.add(lblDistanceBetweenSubplots);
		
		lblSubplotsPerCluster = new JLabel("Sub-plots per Cluster Plot");
		lblSubplotsPerCluster.setEnabled(false);
		lblSubplotsPerCluster.setBounds(221, 493, 200, 50);
		contentPane.add(lblSubplotsPerCluster);
		
		lblStart_Y = new JLabel("Starting coordinate y:");
		lblStart_Y.setEnabled(false);
		lblStart_Y.setBounds(10, 716, 200, 50);
		contentPane.add(lblStart_Y);

		lblNewLabel = new JLabel("Selected File:");
		lblNewLabel.setBounds(10, 47, 82, 30);
		contentPane.add(lblNewLabel);
		
		lbl_SelectColumn = new JLabel("Choose one of the Shapefile Column Names for Strata Selection:");
		lbl_SelectColumn.setBounds(10, 115, 450, 30);
		contentPane.add(lbl_SelectColumn);
		
		lblChooseStrataFor = new JLabel("Choose Strata (one or more) for Sampling:");
		lblChooseStrataFor.setBounds(10, 186, 336, 30);
		contentPane.add(lblChooseStrataFor);
		
		lblSamplingDesign = new JLabel("Sampling Design: ");
		lblSamplingDesign.setBounds(10, 373, 200, 30);
		contentPane.add(lblSamplingDesign);
		
		lblDist_x = new JLabel("Distance between points x:");
		lblDist_x.setEnabled(false);
		lblDist_x.setBounds(10, 493, 200, 50);
		contentPane.add(lblDist_x);
		
		lblStartingPoint = new JLabel("Starting point:");
		lblStartingPoint.setEnabled(false);
		lblStartingPoint.setBounds(10, 602, 200, 50);
		contentPane.add(lblStartingPoint);
		
		lblNewLabel_3 = new JLabel("Cluster Sampling:");
		lblNewLabel_3.setBounds(221, 373, 200, 30);
		contentPane.add(lblNewLabel_3);
		
		lblStart_X = new JLabel("Starting coordinate x:");
		lblStart_X.setEnabled(false);
		lblStart_X.setBounds(10, 663, 200, 50);
		contentPane.add(lblStart_X);
		
		lblDist_y = new JLabel("Distance between points y:");
		lblDist_y.setEnabled(false);
		lblDist_y.setBounds(10, 550, 200, 50);
		contentPane.add(lblDist_y);
		// End Labels
		//-----------------------------------------------------------------------------------------
		
		
		
		
		//-----------------------------------------------------------------------------------------
		// JTable stuff
		// Step 1: model (describe data to be contained in table)
		model = new DefaultTableModel(new String[]{"name", "numPlots"},0); // initialize table with 0 rows
		// (adding rows: this will be done dynamically by clicking on btnAdd)
		// model.addRow(new String[]{"Hamburg", "30"}); // how to add rows
		// Step 2: table
		table = new JTable();
		table.setModel(model);
		JScrollPane scrollPaneTable = new JScrollPane();
		scrollPaneTable.setBounds(235, 227, 326, 105);
		contentPane.add(scrollPaneTable);
		scrollPaneTable.setViewportView(table); // so werden auch die Spaltennamen angezeigt (fallen ohne Darstellung in Scrollpane weg)
		// End JTable Stuff
		//-----------------------------------------------------------------------------------------
		
		
		
		// scrollPaneStrataList: JScrollPane as container for the JList that holds the polygon names
		// initialize outside action method so that element is visible before polygons are ava�lable
		scrollPaneStrataList = new JScrollPane();
		scrollPaneStrataList.setBounds(14, 226, 112, 120);
		contentPane.add(scrollPaneStrataList);
		
		
		
		// btnAdd to move selected items from JList to JTable
		JButton btnAdd = new JButton("Add");
		btnAdd.setBounds(136, 243, 89, 23);
		contentPane.add(btnAdd);
		btnAdd.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				// read selected values from JList
				List<String> selectedValues  = strataNamesList.getSelectedValuesList();
				// write selected values to JTable via model
				for(String value : selectedValues){
					model.addRow(new String[]{value, null}); // rows are added to the underlying TableModel, not the actual JTable !!!
				}
				
			}
		});
		
		
		// btnRemove to remove items from JTable
		JButton btnRemove = new JButton("Remove");
		btnRemove.setBounds(136, 273, 89, 23);
		contentPane.add(btnRemove);
		
		textField_H_numPlotsHorizontal = new JTextField();
		textField_H_numPlotsHorizontal.setEnabled(false);
		textField_H_numPlotsHorizontal.setColumns(10);
		textField_H_numPlotsHorizontal.setBounds(221, 703, 86, 20);
		contentPane.add(textField_H_numPlotsHorizontal);
		
		JLabel lblSubplotsHorizontalLine = new JLabel("Sub-plots per vertical line (only H Clusters):");
		lblSubplotsHorizontalLine.setEnabled(false);
		lblSubplotsHorizontalLine.setBounds(221, 663, 300, 50);
		contentPane.add(lblSubplotsHorizontalLine);
		btnRemove.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				
				int[] selectedRows = table.getSelectedRows();
				
				for(int i = 0; i < selectedRows.length; i++){
					int row = table.getSelectedRow();
					model.removeRow(row);
				}
				
//				for(int i = 0; i < table.getSelectedRows().length; i++){
//					model.removeRow(table.getSelectedRow());
//				}
				
			}
		});
		
	
	}
	
	
	
	private void disableSamplingDesignControls(){
		lblDist_x.setEnabled(false);
		lblDist_y.setEnabled(false);
		lblStartingPoint.setEnabled(false);
		lblStart_X.setEnabled(false);
		lblStart_Y.setEnabled(false);
		textField_DistX.setEnabled(false);
		textField_DistY.setEnabled(false);
		textField_Start_X.setEnabled(false);
		textField_Start_Y.setEnabled(false);
		comboBox_StartingPoint.setEnabled(false);		
	}
	
	private void disableClusterControls(){
		lblClusterDesign.setEnabled(false);
		comboBox_ClusterDesign.setEnabled(false);
		lblSubplotsPerCluster.setEnabled(false);
		textField_NumSubPlotsPerCluster.setEnabled(false);
		lblDistanceBetweenSubplots.setEnabled(false);
		textFieldDistBetweenSubPlots.setEnabled(false);	
		lblSubplotsVerticalLine.setEnabled(false);
		textField_H_numPlotsVertical.setEnabled(false);
	}
}
