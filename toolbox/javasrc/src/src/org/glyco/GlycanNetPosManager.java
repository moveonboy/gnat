package org.glyco;

import java.awt.Rectangle;
import java.awt.geom.Line2D;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.eurocarbdb.application.glycanbuilder.*;

public class GlycanNetPosManager extends Geometry{
	protected GraphicOptions theGraphicOptions;
	protected HashMap<Glycan,Rectangle> glycanPosMap;
	protected HashMap<GlycanRxn,Line2D> glycanRxnPosMap;

	public GlycanNetPosManager()
	{
		this.glycanPosMap = new HashMap<Glycan, Rectangle>();
		this.glycanRxnPosMap = new HashMap<GlycanRxn, Line2D>();
	}

	public GraphicOptions getGraphicOptions() {
		return this.theGraphicOptions;
	}

	public void setGraphicOptions(GraphicOptions opt) {
		this.theGraphicOptions = opt;
	}

	public Iterator<Map.Entry<GlycanRxn, Line2D>> iterator()
	{
		return this.glycanRxnPosMap.entrySet().iterator();
	}

	public Rectangle getGlycanPos(Glycan theGlycan){
		return this.glycanPosMap.get(theGlycan);
	}

	public Line2D getGlycanRxnPos(GlycanRxn theGlyanRxn){
		return this.glycanRxnPosMap.get(theGlyanRxn);
	}

	public void reset() {
		this.glycanPosMap.clear();
		this.glycanRxnPosMap.clear();
	}

	public boolean addGlycanPos(Glycan theGlycan,Rectangle theBox){
		if((theGlycan!=null)&&(theBox!=null)) {
			this.glycanPosMap.put(theGlycan, theBox);
			return true;
		} else{
			return false;
		}
	}

	public boolean addGlycanRxnPos(GlycanRxn theGlycanRxn,Line2D theLine){
		if((theGlycanRxn!=null)&&(theLine!=null))  {
			this.glycanRxnPosMap.put(theGlycanRxn, theLine);
			return true;
		} else{
			return false;
		}
	}

	public void updateRxnsInPathway(GlycanNet glycanNet) {
		//local fields
		double x1,x2,y1,y2;
		for(int i=0;i<glycanNet.getRxns().size();i++){
			Glycan theGlycanReac = glycanNet.getRxns().get(i).getReactants().get(0).getStructure();
			Glycan theGlycanProd = glycanNet.getRxns().get(i).getProducts().get(0).getStructure();
			if((this.glycanPosMap.get(theGlycanReac)==null)|| (this.glycanPosMap.get(theGlycanProd)==null)) continue;

			Glycan leftGlycan,rightGlycan;
			
			if(this.glycanPosMap.get(theGlycanReac).getMinX()<this.glycanPosMap.get(theGlycanProd).getMinX()){
			    leftGlycan  = theGlycanReac;
			    rightGlycan = theGlycanProd;
			}else{
			    leftGlycan = theGlycanProd;
			    rightGlycan = theGlycanReac;
			}			
			
			x1 = this.glycanPosMap.get(leftGlycan).getMaxX();
			y1 = this.glycanPosMap.get(leftGlycan).getCenterY();		    			    
			
			x2 = this.glycanPosMap.get(rightGlycan).getMinX();
			y2 = this.glycanPosMap.get(rightGlycan).getCenterY();
			
			Line2D theLine = new Line2D.Double(x1,y1,x2,y2);
			this.glycanRxnPosMap.put(glycanNet.getRxns().get(i), theLine);
		}
	}
}






