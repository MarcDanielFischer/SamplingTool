
package com.arbonaut.sampling;

import org.opengis.referencing.crs.CoordinateReferenceSystem;
import com.vividsolutions.jts.geom.Geometry;


/**
 * Represents a single sample plot. 
 *  
 * Sample plot stores location/geometry of the plot, coordinate reference system 
 * (crs) and a number of attributes such as the unique plot or cluster numbers, 
 * stratum information etc. 
 *
 * The class is immutable and cannot be modified after created. 
 */
public class SamplePlot {
	
	// properties
	private Geometry geometry;
	private int clusterNr;
	private int plotNr;
	private String stratumName;
	private CoordinateReferenceSystem CRS;     
	private double weight = -1.0;
	

	
	// constructors (3 versions, 2 overloads in order to match variable input param number)
	


	/**
	 * Constructs a Plot object without clusterNr or plotNr.
	 * @param point
	 * @param stratumName
	 * @param CRS
	 */
	public SamplePlot(Geometry geom, String stratumName, CoordinateReferenceSystem CRS) {
		this.geometry = geom;
		this.stratumName = stratumName;
		this.CRS = CRS;
	}
	
	/**
	 * Constructs a Plot object with a plotNr, but without a clusterNr (for CLUSTERSAMPLING_NO option)
	 * @param point
	 * @param stratumName
	 * @param CRS
	 * @param plotNr
	 */
	public SamplePlot(Geometry geom, String stratumName, CoordinateReferenceSystem CRS, int plotNr) {
		this(geom, stratumName, CRS);
		this.plotNr = plotNr;
	}
	
	/**
	 * Constructs a Plot object with plotNr and clusterNr
	 * @param point
	 * @param stratumName
	 * @param CRS
	 * @param plotNr
	 * @param clusterNr
	 */
	public SamplePlot(Geometry geom, String stratumName, CoordinateReferenceSystem CRS, int plotNr, int clusterNr) {
		this(geom, stratumName, CRS, plotNr);
		this.clusterNr = clusterNr;
	}

	
	// Getters 
	public Geometry getGeometry() {
		return geometry;
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

	public CoordinateReferenceSystem getCRS() {
		return CRS;
	}

	public double getWeight() {
		return weight;
	}

}
