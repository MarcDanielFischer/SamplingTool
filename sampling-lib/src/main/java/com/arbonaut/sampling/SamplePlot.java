
package com.arbonaut.sampling;

import org.opengis.referencing.crs.CoordinateReferenceSystem;
import com.vividsolutions.jts.geom.Geometry;


/**
 * Represents a single sample plot. 
 *  
 * Sample plot stores geometry of the plot and a number of attributes such as 
 * the unique plot or cluster numbers, stratum information etc. 
 *
 * The class is immutable and cannot be modified after created. 
 */
public class SamplePlot {
	
    public static final int NO_CLUSTER = -1;
    public static final double NO_WEIGHT = -1;
    
	// properties
    private Geometry geometry;
    private CoordinateReferenceSystem CRS;
    
	private int plotNr = 0;
	private String stratumName = null;
	private int clusterNr = NO_CLUSTER;
	private double weight = NO_WEIGHT;
	

	
	/**
	 * Constructs a SamplePlot object without clusterNr or plotNr.
	 * @param point
	 * @param stratumName
	 * @param CRS
	 */
	public SamplePlot(Geometry geom, CoordinateReferenceSystem crs, String stratumName) {
		this.geometry = geom;
        this.CRS = crs;
		this.stratumName = stratumName;
	}
	
	/**
	 * Constructs a SamplePlot object with a plotNr, but without a clusterNr (for CLUSTERSAMPLING_NO option)
	 * @param point
	 * @param stratumName
	 * @param CRS
	 * @param plotNr
	 */
	public SamplePlot(Geometry geom, CoordinateReferenceSystem crs, String stratumName, int plotNr) {
		this(geom, crs, stratumName);
        this.plotNr = plotNr;
	}
	
	/**
	 * Constructs a SamplePlot object with plotNr and clusterNr
	 * @param point
	 * @param stratumName
	 * @param CRS
	 * @param plotNr
	 * @param clusterNr
	 */
	public SamplePlot(Geometry geom, CoordinateReferenceSystem crs, String stratumName, int plotNr, int clusterNr) {
		this(geom, crs, stratumName, plotNr);
		this.clusterNr = clusterNr;
	}
    
    
    /**
	 * Constructs new SamplePlot from existing instance but with geometry replaced.
	 */
	public SamplePlot(SamplePlot plot, Geometry newgeom, CoordinateReferenceSystem newcrs) {
       this.geometry = newgeom;
       this.CRS = newcrs;
	   this.clusterNr = plot.clusterNr;
	   this.plotNr = plot.plotNr;
	   this.stratumName = plot.stratumName;
	   this.weight = plot.weight;
	}
    

	// Getters 
	public Geometry getGeometry() {
		return geometry;
	}

	public CoordinateReferenceSystem getCRS() {
		return CRS;
	}

	public int getClusterNr() {
		return clusterNr;
	}

	public int getPlotNr() {
		return plotNr;
	}

	public String getStratumName() {
		return stratumName;
	}

	public double getWeight() {
		return weight;
	}

}
