
package com.arbonaut.sampling;

import java.io.File;
import com.arbonaut.sampling.Sampling;
import com.arbonaut.sampling.design.SamplingDesign;
import com.arbonaut.sampling.design.SimpleRandomSample;
import com.arbonaut.sampling.design.plot.PlotGeometry;
import com.arbonaut.sampling.design.plot.PointPlot;
import com.arbonaut.sampling.design.plot.CircularPlot;
import com.arbonaut.sampling.design.cluster.PlotCluster;
import com.arbonaut.sampling.design.cluster.IShapedCluster;
import com.arbonaut.sampling.io.SampleWriter;
import com.arbonaut.sampling.io.CsvSampleWriter;
import com.arbonaut.sampling.io.FilteredSampleWriter;
import com.arbonaut.sampling.io.filter.CenterpointSampleFilter;
import com.arbonaut.sampling.io.filter.ReprojectingSampleFilter;
import org.geotools.referencing.crs.DefaultGeographicCRS;


/*
 * Simple test of the sampling library.
 */ 
public class SamplingApp {
	
    public static void main(String[] args) 
    {
        try
        {        
            // Which plot geometry do we want? Uncomment one option
            //PlotGeometry plotGeom = new PointPlot();
            PlotGeometry plotGeom = new CircularPlot(9);
            
            // What kind of clustering? Uncomment one option
            PlotCluster cluster = null; // No clustering
            //PlotCluster cluster = new IShapedCluster(100, 3);
            
            // Simple random sampling
            SamplingDesign design = new SimpleRandomSample(plotGeom, cluster);      
    
            // We'll write to CSV file        
            SampleWriter csvWriter = new CsvSampleWriter(
                new File("c:\\temp\\sample1.csv"),
                design, 
                false            // If true, writes X,Y columns (requires point geometry)
                );        
                
                   
            // We can apply some filtering before we write to the csv file. 
            // Filters are applied in the order they're added 
                            
            FilteredSampleWriter writer = new FilteredSampleWriter(csvWriter);     

            // Uncomment this reproject the output to WGS84 coordinates
            //writer.addFilter(new ReprojectingSampleFilter(DefaultGeographicCRS.WGS84));    
                            
            // Uncomment this to convert the sample geometry to points
            writer.addFilter(new CenterpointSampleFilter());    

            try
            {
                // Now we could just call design.generate(...) but we'll use the handy
                // methods in Sampling class which do some of the boring work for us  
                
                Sampling.generateSample(
                    new File("c:\\temp\\samplingCensus2.shp"),
                    100,
                    false,   // Force UTM? 
                    design, writer);
            }
            finally
            {
                writer.close();
            }
                
            System.out.println("Done");                
        }
        catch (Exception error)
        {
            System.out.println("Error: " + error.toString());        
            error.printStackTrace();
        }
    }

 
}
