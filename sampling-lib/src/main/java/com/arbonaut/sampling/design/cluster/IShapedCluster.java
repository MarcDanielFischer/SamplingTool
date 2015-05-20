
package com.arbonaut.sampling.design.cluster;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.CoordinateSequence;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;


/**
 * I-shaped plot cluster.
 * I-shaped cluster is a series of plots with an offset applied to the y 
 * coordinate. 
 */
public class IShapedCluster implements PlotCluster {
    
    private double ydistance;
    private int numPlots;  
    
    /**
     * IShapedCluster constructor.
     *
     * @param ydistance  Distance between sub-plots along the y axis     
     * @param numPlots   Number of plots in the cluster     
     */              
    public IShapedCluster(double ydistance, int numPlots)
    {
        this.ydistance = ydistance;
        this.numPlots = numPlots;
    }
    
    public CoordinateSequence create(Coordinate centre)
    {
        Coordinate[] coords = new Coordinate[this.numPlots];
        
        double x = centre.x;
		double y = centre.y;
        
        for (int n = 0; n != this.numPlots; n++)
        {
            coords[n] = new Coordinate(x, y);
            y+= this.ydistance;
        }
        
        return new CoordinateArraySequence(coords);
    
    }
  
}
