package production;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

/**
 * This class provides static methods to write the generated data
 * to output files in different formats.
 * @author daniel
 *
 */
public class OutputUtilities {

	
/**
 * Shows a JFileChooser dialog that prompts the user
 * to specify a file and return that file. To be used e.g. for getting output files.
 * 
 * @return File object
 */
public static File getFile() throws Exception{
	// Problem hier: Methode muss File-Objekt zurückgeben. Wenn User den SaveFileDialog abbricht, ohne eine Datei anzugeben,
	// muss ich hier nach der Programmierlogik eine Exception werfen
	// TODO automatically add ".csv" file extension
	JFileChooser saveFileDialog = new JFileChooser(); // JFileChooser kommt aus javax.swing und hat nichts mit GeoTools zu tun 
	saveFileDialog.setCurrentDirectory(new File("C:\\"));
	//    FileNameExtensionFilter filter = new FileNameExtensionFilter(
	//        "CSV files", "csv");
	//    chooser.setFileFilter(filter);
	//saveFileDialog.setDialogType(JFileChooser.SAVE_DIALOG);
	if(saveFileDialog.showSaveDialog(null) != JFileChooser.APPROVE_OPTION) { // return from method if users cancels the dialog (ie, user does not click on "Save"-Button)
		throw new Exception("SaveFileDialog cancelled without specifying an output file");
	}
	// create output file 
	File saveFile = saveFileDialog.getSelectedFile();
	if(saveFile.exists()){
		int fileOverwriteChoice = JOptionPane.showConfirmDialog(null, "File already exists. Do you want to overwrite the existing file?");
		if(fileOverwriteChoice == JOptionPane.YES_OPTION){
			saveFile.createNewFile();
		}else{
			throw new Exception("you have chosen not to overwrite an existing file");
		}

	}else{
		saveFile.createNewFile();
	}
	return saveFile;
}
	
	/**
	 * Writes the outpot plots generated during the sampling process to a CSV output file.
	 * @param file output file
	 * @param samplePlots output plots to be written to the file
	 * @param clusterSampling specifies whether cluster numbers are to be added to each plot
	 * @param weightedSampling specifies whether weight values are to be added to each plot
	 * @throws Exception
	 */
	public static void writeCSVoutput(File file, ArrayList<Plot> samplePlots, int clusterSampling, boolean weightedSampling) throws Exception{
		// TODO vl doch BufferedWriter nehmen --> ist der irgendwie sicherer?

		FileWriter fileWriter = new FileWriter(file, false); // second param is for appending to file (Yes/No) --> 
		
		// write file header
		// adapt file header for 4 different cluster sampling and weighted sampling option combinations:
		
		// false false
		if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_NO && weightedSampling == false){ 
			fileWriter.write("\"Plot_nr\",\"X\",\"Y\",\"Stratum\"\n");
		}
		// true false
		if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES && weightedSampling == false){ 
			fileWriter.write("\"Cluster_nr\",\"Plot_nr\",\"X\",\"Y\",\"Stratum\"\n");
		}
		// false true
		if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_NO && weightedSampling == true){ 
			fileWriter.write("\"Plot_nr\",\"X\",\"Y\",\"Plot_weight\",\"Stratum\"\n");
		}
		// true true 
		if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES && weightedSampling == true){ 
			fileWriter.write("\"Cluster_nr\",\"Plot_nr\",\"X\",\"Y\",\"Plot_weight\",\"Stratum\"\n");
		}
		
		// write Plot number in the beginning
		
		
		for(Plot plot : samplePlots){
			// write clusterNr only to output if Cluster Sampling is the chosen sampling option
			if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES){
				fileWriter.write(Integer.toString(plot.getClusterNr())+ ",");
			}
			fileWriter.write(Integer.toString(plot.getPlotNr()) + ",");
			double x = plot.getPoint().getCoordinate().x;
			double y = plot.getPoint().getCoordinate().y;
			fileWriter.write(Double.toString(x) + ",");
			fileWriter.write(Double.toString(y) + ",");
			
			if(weightedSampling == true){
				fileWriter.write(Double.toString(plot.getWeight()) + ",");
			}

			fileWriter.write(plot.getStratumName() + "\n");


		}

		fileWriter.close();


		//		// the following two lines are needed to derive appropriate index numbers for the plots
		//		String plotStratumName = null;
		//		int j = 0;
		//		
		//		// iterate over ArrayList and write Point values to output file
		//		for(int i = 0; i < samplePlots.size(); i++){
		//			Point currentPlot = samplePlots.get(i);
		//			double x = currentPlot.getCoordinate().x;
		//			double y = currentPlot.getCoordinate().y;
		//			fileWriter.write(Double.toString(x) + ",");
		//			fileWriter.write(Double.toString(y) + ",");
		//			// the following if-else block is to derive the plot index number for the output file (index starts at one for each stratum)
		//			if (plotStratumName == (String)currentPlot.getUserData()){ // if plotStratumName is same as for the plot before, then continue increasing index j
		//				j++;				
		//			}else{
		//				plotStratumName = (String)currentPlot.getUserData(); // if plotStratumName has changed, then reset index j to 1
		//				j = 1;
		//			}
		//			fileWriter.write((j) + ","); // "Plot_nr" derived from samplePlots-ArrayList index position 
		//			fileWriter.write((String)currentPlot.getUserData() + "\n"); // der Polygonname fehlt jetzt noch
		//		}
		//		fileWriter.close();
	}
	
	/**
	 * TODO to be implemented yet.
	 */
	public static void writeSHPoutput(){
		
	}
	
	/**
	 * TODO to be implemented yet.
	 */
	public static void writeXMLoutput(){
		
	}

}
