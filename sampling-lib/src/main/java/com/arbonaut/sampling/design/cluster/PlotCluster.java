
package com.arbonaut.sampling.design.cluster;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.CoordinateSequence;


/**
 * Plot cluster shape.
 */
public interface PlotCluster {
	
    /**
     * Generates new cluster of plots from some origin point. 
     *       
     * @returns Sequence of coordinates such that each XY pair defines the 
     *          centre of a sample plot in the cluster.     
     */         
    public CoordinateSequence create(Coordinate centre);
  
}
