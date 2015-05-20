
package com.arbonaut.sampling.io;

import java.lang.Exception;
import com.arbonaut.sampling.SamplePlot;

/**
 * Handles persistence of a sample.
 * 
 * SampleWriter implementations provide functionality to write sample to some 
 * persistent format as it is being generated.      
 */
public interface SampleWriter extends java.io.Closeable {

    /**
     * Writes given sample plot to the underlying storage.
     */         
    public void write(SamplePlot plot) throws Exception; 	
	

}
