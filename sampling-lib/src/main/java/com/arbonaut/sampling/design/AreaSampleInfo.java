
package com.arbonaut.sampling.design;

import java.util.Set;
import java.util.Map;
import java.util.HashMap;
import java.util.Collections;
                                                       

/**
 * Contains information about an area sample that can be used to infer its 
 * statistical properties.  
 */
public class AreaSampleInfo implements SampleInfo {

    private class StratumData
    {
        int numPlots;
        double area;
        public StratumData() { numPlots = 0; area = 0; }
    }


    int numPlotsTotal;
    Map<String, StratumData> strata;  
    

    public int dimension() { return 2; }
    public int getNumPlotsTotal() { return numPlotsTotal; }
    public void setNumPlotsTotal(int x) { numPlotsTotal = x; }
    
    public Set<String> getStrataNames() 
    {        
        return Collections.unmodifiableSet(strata.keySet());
    }
    
    public void addStratum(String name)
    {
        if (strata == null) strata = new HashMap<String, StratumData>();
        strata.put(name, new StratumData());
    }
    
    public int getNumPlots(String stratumName) throws IllegalArgumentException
    {
        if (strata == null || !strata.containsKey(stratumName))
            throw new IllegalArgumentException("Stratum "+stratumName+" does not exist");
        return strata.get(stratumName).numPlots; 
    }
    
    public void setNumPlots(String stratumName, int n) throws IllegalArgumentException
    {
        if (strata == null || !strata.containsKey(stratumName))
            throw new IllegalArgumentException("Stratum "+stratumName+" does not exist");
        strata.get(stratumName).numPlots = n;
    }
    
    /**
     * Returns the area of a stratum.
     * @throws IllegalArgumentException if such stratum does not exist      
     */         
    public double getArea(String stratumName) throws IllegalArgumentException
    {
        if (strata == null || !strata.containsKey(stratumName))
            throw new IllegalArgumentException("Stratum "+stratumName+" does not exist");
        return strata.get(stratumName).area; 
    }
    
    /**
     * Sets the area of a stratum.
     * @throws IllegalArgumentException if such stratum does not exist      
     */         
    public void setArea(String stratumName, double area) throws IllegalArgumentException
    {
        if (strata == null || !strata.containsKey(stratumName))
            throw new IllegalArgumentException("Stratum "+stratumName+" does not exist");
        strata.get(stratumName).area = area;
    }
	
    
    
}
