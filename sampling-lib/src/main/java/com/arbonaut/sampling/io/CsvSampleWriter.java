
package com.arbonaut.sampling.io;

import java.io.IOException;
import java.io.File;
import java.io.FileWriter;
import java.util.Collection;
import java.util.ArrayList;

import com.arbonaut.sampling.SamplePlot;
import com.arbonaut.sampling.design.SamplingDesign;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.io.WKTWriter; 


/**
 * Writes sample plots to a csv file.
 * 
 * The writer can write geometry either in WKT format or as X and Y columns (
 * this requires plots with Point geometry).    
 */
public class CsvSampleWriter implements SampleWriter {

    private FileWriter fileWriter = null;
    private WKTWriter wktWriter;

    private String separator;
    private SamplingDesign design;

    // Should only the plot centre coordinate be written as X,Y columns?
    private boolean writeXYColumns;
    
    
     
    /**
     * CsvSampleWriter constructor.
     * 
     * @param file    Output file
     * @param design  Sample design under which the written plots are generated
     * @param writeXYColumns If true, the CSV file will have 'X' and 'Y' columns
     *                       with the plot coordinates. Otherwise 'Geom' column
     *                       will be written with full geometry description in 
     *                       WKT format.                                 
     */             
    public CsvSampleWriter(File file, SamplingDesign design, 
                           boolean writeXYColumns) throws IOException
    {
        this.writeXYColumns = writeXYColumns;
        this.separator = ",";
        this.design = design;
        
        this.fileWriter = new FileWriter(file, false);
          
        // Write the CSV header 
        if (writeXYColumns) {
            fileWriter.write("X"); fileWriter.write(separator);
            fileWriter.write("Y"); 
        }
        else 
            fileWriter.write("Geom"); 
        
        fileWriter.write(separator); fileWriter.write("Plot_nr");
         
        if (design.getPlotClustering() != null) {
            fileWriter.write(separator); fileWriter.write("Cluster_nr");
        } 

        fileWriter.write(separator); fileWriter.write("Stratum");
        fileWriter.write("\n");
        
        wktWriter = new WKTWriter();
        wktWriter.setFormatted(false); 
    }
    
    
    public void write(SamplePlot plot) throws Exception 
    {
        Geometry geom = plot.getGeometry();
        
        // Write center coordinates of the plot
        if (this.writeXYColumns)
        {
            if (!(geom instanceof Point))
                throw new Exception("Sample plot must have point geometry when writing X,Y attributes to csv"); 
        
            Coordinate coord;
            coord = ((Point)geom).getCoordinate();
            
            fileWriter.write(Double.toString(coord.x)); 
            fileWriter.write(separator); 
            fileWriter.write(Double.toString(coord.y));
        }
        // Or the entire geometry as WKT
        else
        {
            writeStr(wktWriter.write(geom), fileWriter);
        }
        
        // Plot number
        fileWriter.write(separator);
        fileWriter.write(Double.toString(plot.getPlotNr()));
        
        // Cluster number
        if (design.getPlotClustering() != null) 
        {
            fileWriter.write(separator);
            fileWriter.write(Double.toString(plot.getClusterNr()));
        }
        
        // Stratum
        fileWriter.write(separator);
        String stratum = plot.getStratumName();
        if (stratum == null || stratum.isEmpty())
            fileWriter.write("\"\"");
        else
            writeStr(stratum, fileWriter);
            
        fileWriter.write("\n");
        
    }
    
    public void close() throws IOException 
    {
        fileWriter.close();
    }	
    
    
    /*
     * Writes string to writer. If string contains spaces, writes it quoted. 
     */     
    private void writeStr(String x, java.io.Writer writer) throws IOException
    {
        boolean quote = x.contains(" "); 
        if (quote) writer.write("\"");
        writer.write(x);
        if (quote) writer.write("\"");
    }
    

}
