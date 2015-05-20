
package com.arbonaut.sampling.util;

import com.arbonaut.sampling.design.plot.PlotGeometry;

import org.opengis.referencing.crs.CoordinateReferenceSystem;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Polygonal;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.prep.PreparedPolygon;
import com.vividsolutions.jts.geom.Location;
import com.vividsolutions.jts.algorithm.locate.PointOnGeometryLocator; 


/**
 * Implements checks whether a sample plot is within a census geometry.  
 * 
 * This is optimized for performance in case of many plots being tested against 
 * a static reference geometry (census/stratum geometry). The way this operates
 * is tied to a particular plot geometry and therefore the same instance should
 * not be used on sample plots whose geometry is created via a different 
 * instance of PlotGeometry.     
 */
public class PlotInPolygonChecker {
	
    PreparedPolygon refGeom;
    PreparedPolygon refGeomShrinked;
    PreparedPolygon refGeomEnlarged;
    
    double plotMBR;
    
    
    /**
     * PlotInPolygonChecker constructor.
     * 
     * @param refGeom  The reference geometry to test against
     * @param plotGeom Plot geometry corresponding to the geometries that will
     *                 be tested                    
     */         
    public PlotInPolygonChecker(Polygonal refGeom, PlotGeometry plotGeom) 
    {
        this.plotMBR = plotGeom.getMBR();
        
        // We'll use JTS PreparedPolygon since the reference geometry does not 
        // change to speed up the tests significantly.   
        
        this.refGeom = new PreparedPolygon(refGeom);
        
        // On top of this, we'll use a 'reduced' and 'enlarged' version of the 
        // reference geometry if MBR is > 0. This will be used for quick 
        // screening before the full geometry test is made.
        
        this.refGeomShrinked = null; 
        this.refGeomEnlarged = null;
        
        if (this.plotMBR > 0)
        {  
            Geometry shrinked = ((Geometry)refGeom).buffer(-this.plotMBR);
            if (shrinked != null && !shrinked.isEmpty() && shrinked.getDimension() >= 2) 
            {
                this.refGeomShrinked = new PreparedPolygon((Polygonal)shrinked);
            }
                        
            Geometry enlarged = ((Geometry)refGeom).buffer(+this.plotMBR);           
            if (enlarged != null && !enlarged.isEmpty() && enlarged.getDimension() >= 2) 
            {
                this.refGeomEnlarged = new PreparedPolygon((Polygonal)enlarged);
            }
        }
    }
    
    /**
     * Tests whether the plot must definitely be inside or outside of the 
     * reference geometry or whether full geometry test is needed. 
     * 
     * @returns Location.INTERIOR if the plot is known to be inside, 
     *          Location.EXTERIOR if the plot is known to be outside,
     *          Location.BOUNDARY if further geometry test is needed
     *          
     * @see com.vividsolutions.jts.geom.Location                                 
     */         
    public int centreTest(Coordinate centre)
    {
        // Special case for MBR == 0 (point) in which we can use this test
        // exclusively
        if (this.plotMBR == 0)
        {
            PointOnGeometryLocator ploc = refGeom.getPointLocator(); 
            int location = ploc.locate(centre);
            return (location == Location.INTERIOR)? 
                Location.INTERIOR : Location.EXTERIOR;
        }
    
        // Check against the geometry that has been shrinked by the MBR distance. 
        // If centre is within, that means the entire plot must also be inside 
        // and no further test is necessary
        if (this.refGeomShrinked != null)
        {
            PointOnGeometryLocator ploc = refGeomShrinked.getPointLocator(); 
            if (ploc.locate(centre) == Location.INTERIOR) 
                return Location.INTERIOR;               
        } 

        // Check against the geometry that has been enlarged by the MBR distance. 
        // If centre is outside, that means the entire plot must also be outside 
        // and no further test is necessary
        if (this.refGeomEnlarged != null)
        {
            PointOnGeometryLocator ploc = refGeomEnlarged.getPointLocator(); 
            if (ploc.locate(centre) == Location.EXTERIOR) 
                return Location.EXTERIOR;               
        } 

        // Out of luck, full geometry test is needed
        return Location.BOUNDARY;
    }
    
    
    
    /**
     * Does full geometry test and returns true if passed geometry is within the 
     * reference geometry.      
     */         
    public boolean contains(Geometry g)
    {
        return this.refGeom.contains(g);
    }
    

    

}
