
package com.arbonaut.sampling.design;

import java.util.Set;
                                                       

/**
 * Contains information about a sample that can be used to infer its statistical
 * properties. 
 * 
 * The info structure is common regardless of the sample design but 
 * differentiated for different sampling dimensions. This abstract base has the
 * common attributes.      
 */
public interface SampleInfo {

    /**
     * Dimension of the sample. Should return either 1 for samples placed along 
     * curves or 2 for area samples. 
     */         
    public int dimension();
    
    /** Returns the total number of plots in all strata. */    
    public int getNumPlotsTotal();

    /** Sets the total number of plots in all strata. */    
    public void setNumPlotsTotal(int x);
    
    /** Returns read-only view of the set of all strata. */         
    public Set<String> getStrataNames();
    
    /**
     * Adds new stratum to the set of strata.
     * Nothing happens if the stratum already exists.     
     */         
    public void addStratum(String name);
    
    /**
     * Returns the number of plots placed in a stratum. 
     * @throws IllegalArgumentException if such stratum does not exist      
     */         
    public int getNumPlots(String stratumName) throws IllegalArgumentException;
    
    /**
     * Sets the number of plots in a stratum.
     * @throws IllegalArgumentException if such stratum does not exist      
     */         
    public void setNumPlots(String stratumName, int n) throws IllegalArgumentException;
        
          
}
