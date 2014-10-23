package org.glyco;
import com.mxgraph.canvas.mxICanvas;
import com.mxgraph.model.mxCell;
import com.mxgraph.swing.mxGraphComponent;
import com.mxgraph.swing.view.mxInteractiveCanvas;
import com.mxgraph.view.mxCellState;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.util.Vector;
import javax.swing.CellRendererPane;
import org.eurocarbdb.application.glycanbuilder.*;


public class GlycanSwingCanvas extends mxInteractiveCanvas {
	protected CellRendererPane rendererPane = new CellRendererPane();
	protected GlycanRenderer theStructureRenderer;
	protected PositionManager thePosManager = new PositionManager();
	protected BBoxManager theBBoxManager = new BBoxManager();
	protected mxGraphComponent graphComponent;
	protected GraphicOptions theGraphicOptions;
	protected static boolean PRESERVE_IMAGE_ASPECT=false;
	public boolean loaded = false;
        public boolean showRedEnd,showLinkage,showMass;

	public GlycanSwingCanvas(mxGraphComponent graphComponent)
	{
		BuilderWorkspace theBuilderWorkspace = new BuilderWorkspace();
		this.theStructureRenderer = theBuilderWorkspace.getGlycanRenderer();
		this.theGraphicOptions = this.theStructureRenderer.getGraphicOptions() ;//theBuilderWorkspace.getGraphicOptions();
		this.theGraphicOptions.MARGIN_LEFT=0;
		this.theGraphicOptions.MARGIN_RIGHT=0;		

		this.graphComponent = graphComponent;
		this.graphComponent.setBackground(Color.white);
	//	vertexRenderer.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
	//	vertexRenderer.setHorizontalAlignment(JLabel.CENTER);
	//	vertexRenderer.setBackground(graphComponent.getBackground().darker());
	//	vertexRenderer.setOpaque(true);				
	}
	
	public BufferedImage loadImage(String image){
		BufferedImage result = super.loadImage(image);
		System.out.println(image+"= "+result);
		return result;		
	}
	
	public GraphicOptions getGraphicOptions(){
		return theGraphicOptions;
	}
	
	public void setGraphicOptions(GraphicOptions theGraphicOptions){
		this.theGraphicOptions = theGraphicOptions;
	}
	
	public double validateGlycanScale(double scale){
		int minSize=10;
		// System.out.println("scale factor for the glycan = "+scale);
		double scaleEff = scale;
		if(this.theGraphicOptions.NODE_SIZE * scale < minSize){
			scaleEff = (double) minSize/(double) this.theGraphicOptions.NODE_SIZE;
			return scaleEff;
		}else{
			return scale;
		}		
	}
	
	public void setGlycanScale(double scale){    
                this.theGraphicOptions.setScale(scale);     
                if(scale<1){
			this.theGraphicOptions.SHOW_REDEND=false;
			this.theGraphicOptions.SHOW_INFO  =false;
                        this.theGraphicOptions.SHOW_MASSES=false;
                        this.theGraphicOptions.setDisplay("compact");
	        }else if(scale<3){
                    if(this.showLinkage){
                        this.theGraphicOptions.setDisplay("linkage");
                    }
                    this.theGraphicOptions.SHOW_REDEND = this.showRedEnd;
                    this.theGraphicOptions.SHOW_MASSES = this.showMass;                        
                }else {
                    this.theGraphicOptions.setDisplay("normalinfo");
                    this.theGraphicOptions.SHOW_REDEND = true;
                    this.theGraphicOptions.SHOW_MASSES = true;
                    this.theGraphicOptions.SHOW_INFO = true;

                }
                this.theGraphicOptions.setScale(scale);     
                       
        }
		
	public void drawVertex(mxICanvas canvas,mxCellState state, String label)
	{

		try{	
			mxCell themxCell = (mxCell) ((mxCell) state.getCell()).getValue();
			GlycanSpecies s = (GlycanSpecies) themxCell.getValue();

			//System.out.println("The Glycan is "+s.toString());		
			Vector<Glycan> a = new Vector<Glycan>();
			a.add(s.getStructure());
			Graphics2D cg;			
			
			boolean show_masses = this.theStructureRenderer.getGraphicOptions().SHOW_MASSES;
			boolean show_redend = this.theStructureRenderer.getGraphicOptions().SHOW_REDEND;
                        //System.out.println("show mass "+show_masses);
                        //System.out.println("show redend "+show_redend);
                        
			theStructureRenderer.getGraphicOptions().MARGIN_TOP=0;
                        theStructureRenderer.getGraphicOptions().MARGIN_LEFT=5;
			theStructureRenderer.getGraphicOptions().MARGIN_BOTTOM=0;
                        theStructureRenderer.getGraphicOptions().MARGIN_RIGHT=5;

			Rectangle bbbox = theStructureRenderer.computeBoundingBoxes(a, 
					show_masses, show_redend, this.thePosManager, this.theBBoxManager,true);
			//state.setRect(bbbox.getX(), bbbox.getY(),bbbox.getWidth(), bbbox.getHeight());

			int x = (int) (state.getX() + translate.x);
			int y = (int) (state.getY() + translate.y);
			int w = (int) (state.getWidth());
			int h = (int) (state.getHeight());
			cg = (Graphics2D) g.create(x, y, w, h);	
			
			this.setGlycanScale(scale);;
                        //this.theStructureRenderer.getGraphicOptions().setScale(this.theStructureRenderer.getGraphicOptions().SCALE);
                        
                        //System.out.println("scale is"+this.theStructureRenderer.getGraphicOptions().SCALE);
			//if(this.theStructureRenderer.getGraphicOptions().SCALE>1) this.theGraphicOptions.SHOW_INFO  =true;
                        
			try {	
				//System.out.println("debugging point 1: "+bbbox.getHeight()+" "+bbbox.getWidth());
				for(int i=0;i<a.size();i++) theStructureRenderer.paint(cg, a.get(i), null,
						null, show_masses, show_redend, this.thePosManager, this.theBBoxManager);   
			} catch(Exception e){
				throw e;
			}
			finally {
				cg.dispose();
			}
		}catch(Exception e){
			return;
		//	super.drawVertex(canvas,state, label);
		//	System.out.println(e.getMessage());
		//	vertexRenderer.setText(label);
		//	rendererPane.paintComponent(g, vertexRenderer, graphComponent,
		//			(int) state.getX() + translate.x, (int) state.getY()
		//			+ translate.y, (int) state.getWidth(),
		//			(int) state.getHeight(), true);
		}
	}

	public GlycanRenderer getGlycanRenderer() {
		// TODO Auto-generated method stub
		return this.theStructureRenderer;
	}

}
