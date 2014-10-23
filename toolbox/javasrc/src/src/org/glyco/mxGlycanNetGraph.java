package org.glyco;
import java.awt.Dimension;
import com.mxgraph.canvas.mxICanvas;
import com.mxgraph.canvas.mxImageCanvas;
import com.mxgraph.model.mxCell;
import com.mxgraph.model.mxIGraphModel;
import com.mxgraph.view.mxCellState;
import com.mxgraph.view.mxGraph;


public class mxGlycanNetGraph extends mxGraph {

	public static int counter=0;
	
	public mxGlycanNetGraph(){
		super();
	}
	
	public mxGlycanNetGraph(mxIGraphModel parammxIGraphModel){
		 super(parammxIGraphModel, null);	
	}
	
	public void fit(Dimension frameDimension){
		double x = frameDimension.getHeight();
		double y = frameDimension.getWidth();
	}
	
    public void drawState(mxICanvas canvas, mxCellState state,
				boolean drawLabel)
		{   
    		String label = (drawLabel) ? state.getLabel() : "";
    		//System.out.println("the label is "+label);
				
			// Indirection for wrapped swing canvas inside image canvas (used for creating
			// the preview image when cells are dragged)
    		boolean isStateCellVertex = getModel().isVertex(state.getCell());
    		boolean isStateCellHasValue = ((mxCell) state.getCell()).getValue() !=null;
    		boolean iscanvasLocal = canvas instanceof GlycanSwingCanvas;
    		
    		
			if (isStateCellVertex 
					&& canvas instanceof mxImageCanvas
					&& ((mxImageCanvas) canvas).getGraphicsCanvas() instanceof GlycanSwingCanvas)
			{
				((GlycanSwingCanvas) ((mxImageCanvas) canvas).getGraphicsCanvas()).drawVertex(canvas,state, label);				
			}
			else if (isStateCellVertex  && isStateCellHasValue && iscanvasLocal)
			{
				boolean isStateCellGlycan=false;boolean isStateCellString=false;
				try{
				  isStateCellGlycan = ((mxCell) ((mxCell) state.getCell()).getValue()).getValue() instanceof org.glyco.GlycanSpecies;
				}catch(Exception e){
				  isStateCellGlycan = false;
				}
				finally{
				  isStateCellString = ((mxCell) state.getCell()).getValue() instanceof java.lang.String;	
				}
				
				if(isStateCellGlycan){
					((GlycanSwingCanvas) canvas).drawVertex(canvas, state, label);
			    }else if(isStateCellString){
			    	super.drawState(canvas,state,drawLabel);
			    }else{
			    	super.drawState(canvas,state,drawLabel);
			    }			    	
			}else {
				//if(canvas==null) System.out.println("canvas case");
				//if(state==null) System.out.println("state case");
				try{
				super.drawState(canvas,state,drawLabel);
				}catch (Exception e){
					return;
				}
				//System.out.println("draw label is "+drawLabel);
				
			}
		}
}
