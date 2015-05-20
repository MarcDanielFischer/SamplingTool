
package com.arbonaut.sampling.io;

import java.io.IOException;
import java.util.List;
import java.util.ArrayList;

import com.arbonaut.sampling.SamplePlot;
import com.arbonaut.sampling.io.filter.SampleFilter;


/**
 * SampleWriter implementation that can appply one or more filters and then 
 * delegate work to another writer. Filters are applied in the order they are
 * added.   
 */
public class FilteredSampleWriter implements SampleWriter 
{
    private SampleWriter writer; // Underlying writer
    private List<SampleFilter> filters; 
    
    public FilteredSampleWriter(SampleWriter writer) 
    {
        this.writer = writer;
        filters = new ArrayList<SampleFilter>();
    }
    
    public void addFilter(SampleFilter filter)
    {
        filters.add(filter);
    } 

    public void write(SamplePlot plot) throws Exception
    {
        for (SampleFilter f : filters)
            plot = f.filter(plot);
    
        this.writer.write(plot);
    } 	

    public void close() throws IOException 
    { 
        writer.close(); 
    }
	

}
