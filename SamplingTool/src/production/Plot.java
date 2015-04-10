/**
 * This class is used to store the output plots created in the sampling process.
 * The main point in creating a class on its own for the output plots is that
 * additional attributes besides the (Point) Geometry can be stored this way.
 * If we would just use the class com.vividsolutions.jts.geom.Point for the plots,
 * Coordinate Reference System information could be stored in the "srid" attribute
 * and even the associated stratum name (the  stratum the plot lies in) could be stuffed in 
 * the "userData" property, but things would get a little tough when it comes to additional attributes
 * like plotNr and clusterNr.  
 * 
 * 
 * The reason we need the CRS property is that the plot must be be reprojected to LatLon before writing it to the output file.
 * The stratum property is needed because we want to write it as additional information to the output file, too.
 */

package pkg_1;

import org.opengis.referencing.crs.CoordinateReferenceSystem;

import com.vividsolutions.jts.geom.Point;

public class Plot {
	
	// properties
	private Point point;
	private int clusterNr;
	private int plotNr;
	private String stratumName;
	private CoordinateReferenceSystem CRS;
	

	
	// constructors (3 versions, 2 overloads in order to match variable input param number)
	
	/**
	 * Constructs a Plot object without clusterNr or plotNr
	 * @param point
	 * @param stratumName
	 * @param CRS
	 */
	public Plot(Point point, String stratumName, CoordinateReferenceSystem CRS) {
		this.point = point;
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
	public Plot(Point point, String stratumName, CoordinateReferenceSystem CRS, int plotNr) {
		this(point, stratumName, CRS);
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
	public Plot(Point point, String stratumName, CoordinateReferenceSystem CRS, int plotNr, int clusterNr) {
		this(point, stratumName, CRS, plotNr);
		this.clusterNr = clusterNr;
	}

	
	
	
	// Getters and Setters
	public Point getPoint() {
		return point;
	}

	public void setPoint(Point point) {
		this.point = point;
	}
	
	public int getClusterNr() {
		return clusterNr;
	}

	public void setClusterNr(int clusterNr) {
		this.clusterNr = clusterNr;
	}

	public int getPlotNr() {
		return plotNr;
	}

	public void setPlotNr(int plotNr) {
		this.plotNr = plotNr;
	}

	public String getStratumName() {
		return stratumName;
	}

	public void setStratumName(String stratumName) {
		this.stratumName = stratumName;
	}

	public CoordinateReferenceSystem getCRS() {
		return CRS;
	}

	public void setCRS(CoordinateReferenceSystem CRS) {
		this.CRS = CRS;
	}


}
