
package production;

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
import javax.swing.JRadioButton;
import java.awt.Font;


@SuppressWarnings("serial")
public class GUI_Designer extends JFrame {
	
	// GUI components
	private JPanel contentPane;
	private JFileChooser openShapeFileDialog;
	private JFileChooser openRasterFileDialog;
	private JButton btnBrowse;
	private JButton btnRunSampling; 
	private JButton btnClose;
	private JButton btnSelectGeotiffRaster;
	private JRadioButton rdbtnWeightedSampling;
	private JTextField textField_Start_X;
	private JTextField textField_Start_Y;
	private JTextField textField_NumSubPlotsPerCluster;
	private JTextField textFieldDistBetweenSubPlots;
	private JTextField textField_H_numPlotsVertical;
	private JTextField textField_DistX;
	private JTextField textField_DistY;
	private JTextField textField_H_numPlotsHorizontal;
	private JTextField textField_BufferSize;
	private JLabel lblSelectedFile;
	private JLabel lblSubplotsVerticalLine;
	private JLabel lblSubplotsHorizontalLine;
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
	private JLabel lblBufferSize;
	private JLabel lblRasterFile;
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
	// alle nicht in jedem Fall benötigten Parameter habe ich schon hier initialisiert,
	// weil (momentan) immer alle Params an die Methode runSampling() übergeben werden
	// und sie dafür initialisiert sein müssen
	File inputShapeFile;
	File inputRasterFile;
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
	double startX = 0; // only for option STARTING_POINT_SPECIFIED 
	double startY = 0; // only for option STARTING_POINT_SPECIFIED 
	int clusterSampling;
	static final int CLUSTER_SAMPLING_NO = 1; // constants to specify clusterSampling
	static final int CLUSTER_SAMPLING_YES = 2;
	int clusterShape = 0; // only for clusterSampling == YES
	static final int I_SHAPE = 1; // constants to specify clusterShape
	static final int L_SHAPE = 2;
	static final int L_SHAPE_UPSIDE_DOWN = 3;
	static final int H_SHAPE = 4;
	static final int SQUARE_SHAPE = 5;
	static final int SQUARE_SHAPE_ROTATED = 6;
	int numSubPlotsinHVerticalLine = 0; // only for clusterShape == H_SHAPE
	int numSubPlotsinHhorizontalLine = 0; // only for clusterShape == H_SHAPE
	int numClusterSubPlots = 0; // only for clusterSampling == YES
	int distBetweenSubPlots = 0; // only for clusterSampling == YES
	double bufferSize;
	boolean weightedSampling = false; // make weighted Sampling not-default 
	
	
	
	
	
	
	

	/**
	 * Launch the application.
	 */
	//TODO evtl in eigene Klassa auslagern, die nur die Anwendung laufen lässt und sonst nichts.
	public static void main(String[] args) {
		// dieses ganze Gedöns mit EventQueue, invokeLater() und run() ist für die Programmausführung durch Threads gedacht
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
		setBounds(100, 100, 684, 902);
		//TODO make GUI window scrollable 
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		// Cancel Button
		btnClose = new JButton("Close");
		btnClose.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				System.exit(0); // TODO evtl. checken, ob das die Option der Wahl ist
			}
		});
		btnClose.setBounds(496, 800, 122, 23);
		contentPane.add(btnClose);
		
		
		
		// "Browse" button
		btnBrowse = new JButton("Select Shapefile");
		btnBrowse.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				openShapeFileDialog = new JFileChooser();
				// TODO change File Dialog start directory to c: oder d:
				openShapeFileDialog.setCurrentDirectory(new File("D:\\_HCU\\_Masterarbeit\\_TestData"));
				if(openShapeFileDialog.showOpenDialog(null) == JFileChooser.APPROVE_OPTION){
					inputShapeFile = openShapeFileDialog.getSelectedFile();
					// write input file path to Label (just for visualization purposes)
					lblSelectedFile.setText(inputShapeFile.toString());
					
					// write input SHP column names to comboBox_StratField immediately after File is selected
					try{
						// clear comboBox_StratField entries first if selected file has been changed
						comboBox_Columns.removeAllItems();
						String[] colNames = SamplingFunctionalityMethods.getSHPColNames(inputShapeFile);
						// add SHP column names to comboBox_StratField dropdown menu
						for(String colName : colNames){
							// TODO sort items alphabetically when adding them to comboBox_StratField
							comboBox_Columns.addItem(colName);
						}
					}catch(Exception e){
						JOptionPane.showMessageDialog(null, "File error: make sure you have selected a Shapefile\n" + e.toString());
					}
				}
			}
		});
		btnBrowse.setBounds(10, 20, 152, 23);
		contentPane.add(btnBrowse);
		
		
		
		// Combo Box "Columns"
		comboBox_Columns = new JComboBox();
		comboBox_Columns.setBounds(10, 100, 200, 30);
		contentPane.add(comboBox_Columns);
		// when a column is selected, all the values in this column will be read automatically and filled as 
        // items into the "Strata" JList
		comboBox_Columns.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				// String actionCommand = e.getActionCommand(); // actionCommand "comboBoxChanged" kommt sowohl bei addItem() als auch bei Item-Auswahl
				JComboBox eventSource = (JComboBox)e.getSource();
				String selectedItem = (String) eventSource.getSelectedItem();
				// Workaround: beim Befüllen der ComboBox wird sofort ein ActionEvent ausgelöst, obwohl das erst passieren soll, wenn der User ein Item selected
				//  --> also über if-Block filtern
				if(selectedItem != "the_geom" && selectedItem != null){ // die Werte für Spalte the_geom sollen nicht ausgelesen werden
					try{
						// read values from specified column and write them to JList
						ArrayList<String> strataNames = SamplingFunctionalityMethods.getColumnValues(inputShapeFile, selectedItem);
						
						strataNamesList = new JList(strataNames.toArray());
						// insert JList into ScrollPane
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
		comboBox_SamplingDesign.addItem(null); // damit die ComboBox am Anfang leer ist --> geht vl noch schöner
		comboBox_SamplingDesign.addItem("Simple Random Sampling");
		comboBox_SamplingDesign.addItem("Systematic Sampling");
		contentPane.add(comboBox_SamplingDesign);
		comboBox_SamplingDesign.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				if((String)comboBox_SamplingDesign.getSelectedItem() == null){ 
					disableSamplingDesignControls();
					
					// delete column "number of Plots" if it has been created before
					if(model.findColumn("Number of Plots") != -1){ // -1 --> column not found
						model.setColumnCount(1); // there is no way of removing a column from a DefaultTableModel directly; by setting the columnCount to 1 all other columns are removed.
					}
				}
				
				if((String)comboBox_SamplingDesign.getSelectedItem() == "Simple Random Sampling"){ 
					disableSamplingDesignControls();
					
					// add column "Number of Plots" only if it does not yet exist
					if(model.findColumn("Number of Plots")== -1){ // -1 --> column not found
					model.addColumn("Number of Plots");
					}
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

					// delete column "Number of Plots" if it has been created before
					if(model.findColumn("Number of Plots") != -1){ // -1 --> column not found
						model.setColumnCount(1); // there is no way of removing a column from a DefaultTableModel directly; by setting the columnCount to 1 all other columns are removed.
					}
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
					// --> die folgenden 4 Zeilen könnte ich evtl. in eine MEthode auslagern, weil sie doppelt vorkommen
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
		comboBox_ClusterDesign.addItem("L Upside Down");
		comboBox_ClusterDesign.addItem("H");
		comboBox_ClusterDesign.addItem("Square");
		comboBox_ClusterDesign.addItem("Rotated Square");
		contentPane.add(comboBox_ClusterDesign);
		comboBox_ClusterDesign.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				// if Cluster Shape == H then enable "number of Plots in one vertical line" Control else disable said control
				if(comboBox_ClusterDesign.getSelectedItem()=="H"){
					lblSubplotsVerticalLine.setEnabled(true);
					textField_H_numPlotsVertical.setEnabled(true);
					lblSubplotsHorizontalLine.setEnabled(true);
					textField_H_numPlotsHorizontal.setEnabled(true);
					lblSubplotsPerCluster.setEnabled(false);
					textField_NumSubPlotsPerCluster.setEnabled(false);
					
					
				}else{
					lblSubplotsVerticalLine.setEnabled(false);
					textField_H_numPlotsVertical.setEnabled(false);
					lblSubplotsHorizontalLine.setEnabled(false);
					textField_H_numPlotsHorizontal.setEnabled(false);
					lblSubplotsPerCluster.setEnabled(true);
					textField_NumSubPlotsPerCluster.setEnabled(true);
				}
			}
		});
		

		// "Run Sampling" Button 
		btnRunSampling = new JButton("Run Sampling");
		btnRunSampling.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent arg0) {
				// alle eingetragenen Samplingparameter zusammentragen und auf die am Klassenanfang deklarierten Variablen schreiben
				// TODO besseres Verfahren überlegen, die ganzen Samplingparameter zu übergeben (zu viele für eine einzelne Methode) --> Hashmap, eigenes Objekt nur für die Parameter?

				//----------------------------------------
				// Start Samplingparameter
				// --> Erfassung der Parameter enthält Samplinglogik

				// TODO see if this type of bulletproofing the code can be achieved in a better way
				boolean allParamsOK = true; // flag variable
				
				if(inputShapeFile == null){
					JOptionPane.showMessageDialog(null, "Make sure you have selected a Shapefile");
					allParamsOK = false;

				}

				sampleColumn = (String)comboBox_Columns.getSelectedItem();

				numStrata = model.getRowCount();
				if(numStrata == 0 && allParamsOK == true){
					JOptionPane.showMessageDialog(null, "Please select at least one stratum");
					allParamsOK = false;
				}

				if(allParamsOK){
					selectedStrata = new String[numStrata];
					// read stratum names from JTable
					for(int i = 0; i < numStrata; i++){ // iterate over rows containing stratum names
						selectedStrata[i] = (String) model.getValueAt(i, 0); // i = row; 0 = "name" column
					}
				}

				// Sampling Design params
				if((String)comboBox_SamplingDesign.getSelectedItem()== null && allParamsOK ){
					JOptionPane.showMessageDialog(null, "Please choose a Sampling Design Option");
					allParamsOK = false;
				}else if((String)comboBox_SamplingDesign.getSelectedItem()=="Simple Random Sampling"){
					samplingDesign = SIMPLE_RANDOM_SAMPLING;
					// Number of plots already specified in Array "numPlotsToBeSampled"
					numPlotsToBeSampled = new int[numStrata];

					try{
						for(int i = 0; i < numStrata; i++){ // iterate over rows
							numPlotsToBeSampled[i] = Integer.parseInt((String)model.getValueAt(i, 1)); // i = row; 1 = "number of Plots" column
						}
					}catch(Exception e){
						JOptionPane.showMessageDialog(null, "Please make sure to type the number of plots to be sampled for each stratum\n and then press Enter\n" + e.toString());
						allParamsOK = false;
					}

				}else if ((String)comboBox_SamplingDesign.getSelectedItem()=="Systematic Sampling"){
					samplingDesign = SYSTEMATIC_SAMPLING;

					try{
						gridDistX = Integer.parseInt(textField_DistX.getText()); 
						gridDistY = Integer.parseInt(textField_DistY.getText()); 
					}catch(Exception e){
						JOptionPane.showMessageDialog(null, "Please make sure to specify grid distances X and Y\n" + e.toString());
						allParamsOK = false;
					}

					// Starting Point 
					if((String)comboBox_StartingPoint.getSelectedItem() == "Specified"){
						startingPoint = STARTING_POINT_SPECIFIED;

						// make sure the TextField input numbers are separated by a point instead of a comma (numbers using decimal commas cannot be parsed by Double.parseDouble())
						String x = textField_Start_X.getText().replace(",", ".");

						String y = textField_Start_Y.getText().replace(",", ".");

						try{
							startX = Double.parseDouble(x);
							startY = Double.parseDouble(y);
						}catch(Exception e){
							JOptionPane.showMessageDialog(null, "Please make sure to specify starting coordinates  X and Y\n" + e.toString());
							allParamsOK = false;
						}


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
					if((String)comboBox_ClusterDesign.getSelectedItem() == null && allParamsOK){
						JOptionPane.showMessageDialog(null, "Please make sure to select a Cluster Design Option");
						allParamsOK = false;
					}

					if(allParamsOK){
						switch((String)comboBox_ClusterDesign.getSelectedItem()){
						case "I": clusterShape = I_SHAPE;
						break;
						case "L": clusterShape = L_SHAPE;
						break;
						case "L Upside Down": clusterShape = L_SHAPE_UPSIDE_DOWN;
						break;
						case "H": clusterShape = H_SHAPE;
						break;
						case "Square": clusterShape = SQUARE_SHAPE;
						break;
						case "Rotated Square": clusterShape = SQUARE_SHAPE_ROTATED;
						break;
						}
					}

					// Number of Plots in a vertical line of the H-Shape
					if(clusterShape == H_SHAPE && allParamsOK){
						try{
							numSubPlotsinHVerticalLine = Integer.parseInt(textField_H_numPlotsVertical.getText());
							numSubPlotsinHhorizontalLine = Integer.parseInt(textField_H_numPlotsHorizontal.getText());
						}catch(Exception e){
							JOptionPane.showMessageDialog(null, "Please specify the number of Sub-Plots in the vertical and horizontal lines of the H-Clusters");
							allParamsOK = false;
						}



					}
					// number of sub-plots per cluster (not for H-shaped clusters)
					if(clusterShape != H_SHAPE && clusterShape != SQUARE_SHAPE_ROTATED && allParamsOK ){
						try{
							numClusterSubPlots = Integer.parseInt(textField_NumSubPlotsPerCluster.getText());
						}catch(Exception e){
							JOptionPane.showMessageDialog(null, "Please select a number of Sub-plots per Cluster Plot");
							allParamsOK = false;
						}

					}

					// Distance between Cluster SubPlots
					if(allParamsOK ){
						try{
							distBetweenSubPlots = Integer.parseInt(textFieldDistBetweenSubPlots.getText());
						}catch(Exception e){
							JOptionPane.showMessageDialog(null, "Please specify the distance between Sub-plots in a Cluster Plot");
							allParamsOK = false;
						}
					}



				}else clusterSampling = CLUSTER_SAMPLING_NO; // use CLUSTER_SAMPLING_NO as default option, unless CLUSTER_SAMPLING_YES is explicitly specified

				// size of the negative Buffer to be applied to strata polygons
				try{
					bufferSize = Double.parseDouble(textField_BufferSize.getText());
				}catch(Exception e){
					JOptionPane.showMessageDialog(null, "Please specify a Plot Radius\n" + e.toString());
					allParamsOK = false;
				}


				// Ende Samplingparameter
				//----------------------------------------


				// Methode erst hinterher aufrufen mit gesammelten Sachen drin --> Methodensignatur ändern
				// TODO Methode vl mit weniger Parameter hinbekommen --> evtl. Params als eigene Objektklasse
				if(allParamsOK ){ 
					try{
						SamplingFunctionalityMethods.runSampling(inputShapeFile, sampleColumn, selectedStrata, samplingDesign, numPlotsToBeSampled, gridDistX, gridDistY, startingPoint, startX, startY, clusterSampling, clusterShape, numSubPlotsinHVerticalLine, numSubPlotsinHhorizontalLine, numClusterSubPlots, distBetweenSubPlots, bufferSize, weightedSampling, inputRasterFile ); 
					}catch(Exception e){
						JOptionPane.showMessageDialog(null, e.toString());
					}	
				}
			}
		}
				);
		btnRunSampling.setBounds(4, 800, 122, 23);
		contentPane.add(btnRunSampling);


		//-----------------------------------------------------------------------------------------
		// JTable stuff
		// Step 1: model (describe data to be contained in table)
//		model = new DefaultTableModel(new String[]{"name", "numPlots"},0); // initialize table with 0 rows
		model = new DefaultTableModel();
		model.addColumn("Stratum Name");
		// (adding rows: this will be done dynamically by clicking on btnAdd)
		// model.addRow(new String[]{"Hamburg", "30"}); // how to add rows
		// Step 2: table
		table = new JTable();
		table.setModel(model);
		JScrollPane scrollPaneTable = new JScrollPane();
		scrollPaneTable.setBounds(389, 226, 229, 105);
		contentPane.add(scrollPaneTable);
		scrollPaneTable.setViewportView(table); // so werden auch die Spaltennamen angezeigt (fallen ohne Darstellung in Scrollpane weg)
		// End JTable Stuff
		//-----------------------------------------------------------------------------------------
		
		
		
		// scrollPaneStrataList: JScrollPane as container for the JList that holds the polygon names
		// initialize outside action method so that element is visible before polygons are avaílable
		scrollPaneStrataList = new JScrollPane();
		scrollPaneStrataList.setBounds(10, 170, 270, 204);
		contentPane.add(scrollPaneStrataList);
		
		
		
		// btnAdd to move selected items from JList to JTable
		JButton btnAdd = new JButton("Add");
		btnAdd.setBounds(290, 243, 89, 23);
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
		btnRemove.setBounds(290, 273, 89, 23);
		contentPane.add(btnRemove);
		btnRemove.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				
				int[] selectedRows = table.getSelectedRows();
				
				for(int i = 0; i < selectedRows.length; i++){
					int row = table.getSelectedRow();
					model.removeRow(row);
				}
			}
		});
		
		
		//-----------------------------------------------------------------------------------------
				// TextFields
				
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
				
				// TextField "H_numPlotsHorizontal"
				textField_H_numPlotsHorizontal = new JTextField();
				textField_H_numPlotsHorizontal.setEnabled(false);
				textField_H_numPlotsHorizontal.setColumns(10);
				textField_H_numPlotsHorizontal.setBounds(221, 731, 86, 20);
				contentPane.add(textField_H_numPlotsHorizontal);
				
				// TextField "BuffferSize"
				textField_BufferSize = new JTextField();
				textField_BufferSize.setText("9");
				textField_BufferSize.setColumns(10);
				textField_BufferSize.setBounds(449, 431, 86, 20);
				contentPane.add(textField_BufferSize);
				// End TextFields
				//-----------------------------------------------------------------------------------------
				
				
				//-----------------------------------------------------------------------------------------
				// Labels
				lblSelectedFile = new JLabel("No File selected yet");
				lblSelectedFile.setBounds(92, 40, 448, 30);
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
				lblNewLabel.setBounds(10, 40, 82, 30);
				contentPane.add(lblNewLabel);
				
				lbl_SelectColumn = new JLabel("Choose one of the Shapefile Column Names for Strata Selection:");
				lbl_SelectColumn.setBounds(10, 70, 450, 30);
				contentPane.add(lbl_SelectColumn);
				
				lblChooseStrataFor = new JLabel("Choose Strata (one or more) for Sampling:");
				lblChooseStrataFor.setBounds(10, 141, 336, 30);
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
				
				lblSubplotsHorizontalLine = new JLabel("<html><body>Sub-plots per horizontal line (only H Clusters):<br> Note: Plots located both on horizontal AND vertical lines are considered to belong only to vertical lines</body></html>");
				lblSubplotsHorizontalLine.setEnabled(false);
				lblSubplotsHorizontalLine.setBounds(221, 680, 300, 50);
				contentPane.add(lblSubplotsHorizontalLine);
				
				lblBufferSize = new JLabel("<html><body>Plot Radius<br>(minimum Distance <br>to Stratum boundary):</body></html>");
				lblBufferSize.setBounds(449, 376, 229, 53);
				contentPane.add(lblBufferSize);
				
				
////////////////////////////////////////////////////////////////////////////////////////////////////
// neu weighted Sampling
				lblRasterFile = new JLabel("");
				lblRasterFile.setFont(new Font("Tahoma", Font.PLAIN, 10));
				lblRasterFile.setEnabled(false);
				lblRasterFile.setBounds(300, 65, 360, 14);
				contentPane.add(lblRasterFile);
				
				btnSelectGeotiffRaster = new JButton("Select GeoTIFF Raster File");
				btnSelectGeotiffRaster.setEnabled(false);
				btnSelectGeotiffRaster.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						openRasterFileDialog = new JFileChooser();
						// TODO change File Dialog start directory to c: oder d:
						openRasterFileDialog.setCurrentDirectory(new File("D:\\_HCU\\_Masterarbeit\\_TestData"));
						if(openRasterFileDialog.showOpenDialog(null) == JFileChooser.APPROVE_OPTION){
							inputRasterFile = openRasterFileDialog.getSelectedFile();
							// write input file path to Label (just for visualization purposes)
							lblRasterFile.setText(inputRasterFile.toString());
							
						}
						
					}
				});
				btnSelectGeotiffRaster.setBounds(378, 33, 183, 23);
				contentPane.add(btnSelectGeotiffRaster);
				

				
				rdbtnWeightedSampling = new JRadioButton("Weighted Sampling");
				rdbtnWeightedSampling.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if(rdbtnWeightedSampling.isSelected()){
							btnSelectGeotiffRaster.setEnabled(true);
							lblRasterFile.setEnabled(true);
							weightedSampling = true;
						}else{
							btnSelectGeotiffRaster.setEnabled(false);
							lblRasterFile.setEnabled(false);
							weightedSampling = false;
						}
						
					}
				});
				rdbtnWeightedSampling.setBounds(389, 10, 128, 23);
				contentPane.add(rdbtnWeightedSampling);
////////////////////////////////////////////////////////////////////////////////////////////////////				

				// End Labels
				//-----------------------------------------------------------------------------------------
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
