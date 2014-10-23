package org.glyco;
import com.mxgraph.canvas.mxICanvas;
import com.mxgraph.canvas.mxSvgCanvas;
import com.mxgraph.io.graphml.mxGraphMlPort;
import com.mxgraph.layout.hierarchical.mxHierarchicalLayout;
import com.mxgraph.layout.*;
import com.mxgraph.layout.orthogonal.mxOrthogonalLayout;
import com.mxgraph.model.mxCell;
import com.mxgraph.model.mxGraphModel;
import com.mxgraph.swing.handler.mxRubberband;
import com.mxgraph.swing.mxGraphComponent;
import com.mxgraph.swing.mxGraphOutline;
import com.mxgraph.util.mxCellRenderer.CanvasFactory;
import com.mxgraph.util.mxEventSource.mxIEventListener;
import com.mxgraph.util.mxUndoableEdit.mxUndoableChange;
import com.mxgraph.util.*;
import com.mxgraph.view.mxGraph;
import com.mxgraph.view.mxLayoutManager;
import com.mxgraph.view.mxStylesheet;
import com.mxgraph.view.mxSwimlaneManager;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;
import org.eurocarbdb.application.glycanbuilder.GraphicOptions;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.SBMLDocument;
import org.sbml.jsbml.SBMLReader;

public class GlycanNetGraphEditor extends JFrame 
{

	/*
	 * 
	 */
	private static final String title ="Glycosylation Network Visualizer";
	private static final long serialVersionUID = -844106998814982739L;

	JMenu fileMenu,importMenu,exprotMenu,viewMenu,scaleMenu;
	JMenuItem exitMenu;
	JMenuBar menuBar;
	JMenuItem menuItem1, menuItem2, menuItem3, UndoMenuItem;	

	private static final String metalLook = "javax.swing.plaf.metal.MetalLookAndFeel";
	private static final String motifLook = "com.sun.java.swing.plaf.motif.MotifLookAndFeel";
	private static final String GTKLook = "com.sun.java.swing.plaf.gtk.GTKLookAndFeel";
	private static final String WindowLook = "com.sun.java.swing.plaf.windows.WindowsLookAndFeel";
	private static final String prefix_KEGG = "http://www.genome.jp/dbget-bin/www_bget?gl:";
	private static final String prefix_GLYCOMEDB = "http://www.glycome-db.org/database/" +
			"showStructure.action?glycomeId=";
        private static final String link_GLYCOMEDB = "http://www.glycome-db.org/About.action";

	GlycanNetRead glycanNetModel;
	Vector<GlycanNet> theGlycanNets;

//	static{
//		System.loadLibrary("sbmlj");
//	}

	/**
	 * 
	 */
	protected mxUndoManager undoManager= new mxUndoManager();

	/**
	 * 
	 */
	protected mxIEventListener undoHandler = new mxIEventListener()
	{
		public void invoke(Object source, mxEventObject evt)
		{
			undoManager.undoableEditHappened((mxUndoableEdit) evt
					.getProperty("edit"));
		}
	};

	/**
	 * 
	 */
	protected mxIEventListener changeTracker = new mxIEventListener()
	{
		public void invoke(Object source, mxEventObject evt)
		{
			setModified(true);
		}
	};
	private boolean modified;


	/**
	 * 
	 */
	public mxGraphOutline getGraphOutline()
	{
		return graphOutline;
	}

	/**
	 * 
	 * @param modified
	 */
	public void setModified(boolean modified)
	{
		boolean oldValue = this.modified;
		this.modified = modified;

		firePropertyChange("modified", oldValue, modified);	
	}

	HashMap<GlycanSpecies,Object> glycanVertexManager = new HashMap<GlycanSpecies,Object>();
	HashMap<GlycanNet,Object> glycanNetVertexManager = new HashMap<GlycanNet,Object>();
	HashMap<GlycanRxn,Object> glycanRxnEdgeManager = new HashMap<GlycanRxn,Object>();
	mxGraphModel model = new mxGraphModel();
	mxGraph graph = new mxGlycanNetGraph(model);
	mxSwimlaneManager swimlaneManager;//r=new mxSwimlaneManager(graph);
	mxGraphComponent graphComponent;	
	mxHierarchicalLayout netLayout;
	GlycanSwingCanvas glycanCanvas;
	GraphicOptions theGraphicsOptions;
	Dimension preferredSize;
	String swimlaneStyle;
	mxStylesheet graphStyleSheet;
	protected Dimension defaultSize; 
	protected File currentFile;
	protected double scale;
	protected JFileChooser fc;
	private JMenu layoutMenu;
	private JMenuItem circleLayoutItem;
	private JMenuItem hieraLayoutItem;
	private JMenuItem compactTreeLayoutItem;
	private JMenu helpMenu;
	private JMenuItem aboutMenu;
	public mxLayoutManager layoutManager;
	private mxGraphOutline graphOutline;

	/**
	 * 
	 */
	protected mxUndoManager createUndoManager()
	{
		return new mxUndoManager();
	}
        
        public static void createAndShowGUI(final Vector<GlycanNet> thePathway,final String frameDimensionChoice,final Dimension frameDimension,
	    final boolean showRedEnd, final boolean showMass, final boolean showLinkage,final String netLayout,final HashMap<String,Integer> glycomeIDs){
            createAndShowGUI(thePathway,frameDimensionChoice,frameDimension,showRedEnd,showMass,showLinkage,netLayout,glycomeIDs,false);
        }
        
         public static void createAndShowGUI(final Vector<GlycanNet> thePathway,final String frameDimensionChoice,final Dimension frameDimension,
	    final boolean showRedEnd, final boolean showMass, final boolean showLinkage,final String netLayout,final HashMap<String,Integer> glycomeIDs,final boolean isfit){
            
             java.awt.EventQueue.invokeLater(new Runnable() {
			public void run() {                            
                                if(thePathway.size()==0){
                                   try {
                                       throw new Exception("No Pathway Provided");
                                   } catch (Exception ex) {
                                       Logger.getLogger(GlycanNetGraphEditor.class.getName()).log(Level.SEVERE, null, ex);
                                    }                               
                                }                            
                            
				GlycanNetGraphEditor frame = new GlycanNetGraphEditor();
                                
                                try {
                                   frame.glycanNetModel=new GlycanNetRead(thePathway);
                                } catch (Exception ex) {
                                   Logger.getLogger(GlycanNetGraphEditor.class.getName()).log(Level.SEVERE, null, ex);
                                }

				//set up menus
				frame.setMenuBar();
				frame.installListeners();

				//set up look and feel
				frame.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

				// create cell styles
				frame.createEdgeStyle();
				frame.createSwimLaneStyle();

			     //frame.setSize(frameDimension.width, frameDimension.height);
                             //  System.out.println("width ="+frameDimension.width);
                             //   System.out.println("height ="+frameDimension.height);

				try {
					frame.visNetFromPathway(thePathway,frameDimensionChoice,frameDimension,
                                                showRedEnd,showMass,showLinkage,netLayout,glycomeIDs,isfit);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}	

				//String filename = "try1.png";
				//System.out.println("export image to file");
				///frame.export(filename);
			}
		}
             );		
            
        }
         
         public void visNetFromPathway(Vector<GlycanNet> thePathway,String frameDimensionChoice,Dimension frameDimension,
	    boolean showRedEnd, boolean showMass, boolean showLinkage,String netLayout,HashMap<String,Integer> glycomeIDs) throws Exception{
              visNetFromPathway(thePathway,frameDimensionChoice,frameDimension,showRedEnd,showMass,showLinkage,netLayout,glycomeIDs,false);
         }

        public void visNetFromPathway(Vector<GlycanNet> thePathway,String frameDimensionChoice,Dimension frameDimension,
	    boolean showRedEnd, boolean showMass, boolean showLinkage,String netLayout,HashMap<String,Integer> glycomeIDs,boolean isfit) throws Exception{
            
                this.setFrameOptions(frameDimensionChoice, frameDimension, showRedEnd, showMass, showLinkage);		
                this.theGlycanNets = thePathway;
                              
              	this.assignGlycomeID(this.theGlycanNets,glycomeIDs);
		this.visNet(this.theGlycanNets);

		//frame.visNet(frame.theGlycanNets);
                new mxRubberband(this.graphComponent);	
		
		///System.out.println("size of the network is "+glycanNetVertexManager.size());
		//this.layoutGraph(this.graph.getDefaultParent(),this.createMxGraphLayout("CompactTree",this.graph));
                this.layoutGlycanNetGraph(netLayout,isfit);
		//System.out.println("the default size for the network is " + 
                //        this.defaultSize.getHeight() + " " + this.defaultSize.getWidth());
		//Dimension frameDimn = new Dimension(this.getWidth(),this.getHeight());
	//	System.out.println("the frame size for the network is "+frameDimn.height+" "+frameDimn.width);
		//if(isfit) this.fit(frameDimn);
                this.setVisible(true);			                   
        };

    protected void setFrameOptions(String frameDimensionChoice, Dimension frameDimension, boolean showRedEnd, boolean showMass, boolean showLinkage) {
        this.setBackground(Color.white);
        this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);		
        this.swimlaneStyle = this.setSwimLaneStype();
        //	this.graphComponent.setBackground(Color.WHITE);
        //	this.graphComponent.setOpaque(false);
        this.graphComponent.getViewport().setOpaque(true);
        this.graphComponent.getViewport().setBackground(Color.WHITE);
        this.getContentPane().setBackground(Color.WHITE);
        this.getContentPane().add(this.graphComponent);

        //handle input arguments

        //second: frame size
        if(frameDimensionChoice.equalsIgnoreCase("Custom")){
                this.setSize(frameDimension);
        }else if(frameDimensionChoice.equalsIgnoreCase("Detailed")){
                this.setSize(this.computeDefaultSize());
        }else if(frameDimensionChoice.equalsIgnoreCase("Preferred")){
                this.setSize(this.computePreferredSize());
        }else {
                this.setSize(this.computePreferredSize());
        }

        //third graphic options
        this.theGraphicsOptions.SHOW_REDEND=showRedEnd;
        this.theGraphicsOptions.SHOW_REDEND_CANVAS=showRedEnd;

        this.theGraphicsOptions.SHOW_MASSES=showMass;
        this.theGraphicsOptions.SHOW_MASSES_CANVAS=showMass;

        if(showLinkage) 
           this.theGraphicsOptions.setDisplay("Linkage");
        else
           this.theGraphicsOptions.setDisplay("compact");
        
        this.glycanCanvas.showLinkage = showLinkage;
        this.glycanCanvas.showRedEnd  = showRedEnd;
        this.glycanCanvas.showMass    = showMass;
        

        // Adds rubberband selection		
        //System.out.println("ready for visulization");
        
        this.graph.setAutoOrigin(true);
        this.graph.setAutoSizeCells(true);	
        this.graph.setCellsEditable(true);
        this.graph.setCellsResizable(true);
    }
      
    public static void createAndShowGUI(final String SBMLfileName,final String glycanFileFormat,final String frameDimensionChoice, 
                final Dimension frameDimension,final boolean showRedEnd, final boolean showMass, 
                final boolean showLinkage,final String netLayout,final HashMap<String,Integer> glycomeIDs){
        createAndShowGUI(SBMLfileName,glycanFileFormat,frameDimensionChoice,frameDimension,showRedEnd,showMass,showLinkage,netLayout,glycomeIDs,false);
    }
    
    
	/**
	 * 
	 */
	public static void createAndShowGUI(final String SBMLfileName,final String glycanFileFormat,final String frameDimensionChoice, 
                final Dimension frameDimension,final boolean showRedEnd, final boolean showMass, 
                final boolean showLinkage,final String netLayout,final HashMap<String,Integer> glycomeIDs,final boolean isfit){
            
		java.awt.EventQueue.invokeLater(new Runnable() {
			public void run() {
				GlycanNetGraphEditor frame = new GlycanNetGraphEditor();	

				//set up menus
				frame.setMenuBar();
				frame.installListeners();

				//set up look and feel
				frame.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

				// create cell styles
				frame.createEdgeStyle();
				frame.createSwimLaneStyle();

				frame.setSize(frameDimension.height, frameDimension.width);

				try {
					frame.visNetFromSBMLFile(SBMLfileName,glycanFileFormat,frameDimensionChoice,frameDimension, 
							showRedEnd,showMass,showLinkage,netLayout,glycomeIDs,isfit);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}	

				mxGraphMlPort port= new mxGraphMlPort(SBMLfileName);

			//	String filename = "try1.png";
			//	System.out.println("export image to file");
			//	frame.export(filename);
			}
		});		
	}

	/**
	 * 
	 */
	public static void createAndShowGUI(final String SBMLfileName,final String glycanFileFormat,final String frameDimensionChoice, 
         final Dimension frameDimension,final boolean showRedEnd, final boolean showMass, final boolean showLinkage,final String netLayout){	
		createAndShowGUI(SBMLfileName,glycanFileFormat,frameDimensionChoice,frameDimension,showRedEnd,showMass,showLinkage,netLayout,null);                	
	}

	public static void main(String[] args) throws Exception
	{   
	  //      final String SBMLfileName="glycoNetExample.xml";
        //        final String SBMLfileName="gnat_test_wtannot.xml";
              //  final String SBMLfileName="OlinkedModel_wtannot.xml";
                  final String SBMLfileName="illustractionforpaper.xml";  
          //      final String SBMLfileName="ub1997_onecompt.xml";
          //      final String SBMLfileName="ub1997.xml";
		final String glycanFileFormat="glycoct_xml";
		final String frameDimensionChoice="Custom";
		final Dimension frameDimension = new Dimension(600,800); 
		final boolean showRedEnd = true;
		final boolean showMass = false;
		final boolean showLinkage = true;
                final boolean isFrameFit = false;
		final String netLayout = "CompactTree";
		
                boolean isFileInput = true;
                
                GlycanNetRead theGlycanNetSBMLRead = new GlycanNetRead();
                theGlycanNetSBMLRead.readSBMLfromFile(SBMLfileName);
                Vector<GlycanNet> theGlycanNets = theGlycanNetSBMLRead.theGlycanRxnNets;
               
                int numSpecies = getSpeciesNum(SBMLfileName);
                
	        final HashMap<String,Integer> glycomeDBIDs = new HashMap();
		for(int i=0;i<numSpecies;i++){
	           glycomeDBIDs.put(theGlycanNetSBMLRead.theGlycanArray.get(i).toGlycoCT(),-1);   		
		}
                
                if(!isFileInput){
                 createAndShowGUI(theGlycanNets, frameDimensionChoice, frameDimension, 
			showRedEnd, showMass, showLinkage, netLayout, glycomeDBIDs, isFrameFit);           
                }else{
                 createAndShowGUI(SBMLfileName,  glycanFileFormat, frameDimensionChoice, frameDimension, 
			showRedEnd, showMass, showLinkage, netLayout, glycomeDBIDs, isFrameFit);  
                }
                
	}  
        
        public static int getSpeciesNum(String sbmlFileName) throws Exception {
		SBMLReader theSBMLReader = new SBMLReader();
		SBMLDocument thesbmlDocument = theSBMLReader.readSBMLFromFile(sbmlFileName);
		Model theSBMLModel = thesbmlDocument.getModel();
		return (int) theSBMLModel.getNumSpecies();		
	}

	protected void installListeners()
	{
		// Installs mouse wheel listener for zooming
		MouseWheelListener wheelTracker = new MouseWheelListener()
		{
			/**
			 * 
			 */
			public void mouseWheelMoved(MouseWheelEvent e)
			{
				if (e.getSource() instanceof mxGraphOutline
						|| e.isControlDown())
				{
					GlycanNetGraphEditor.this.mouseWheelMoved(e);

				}
			}
		};


		graphOutline.addMouseWheelListener(wheelTracker);
		graphComponent.addMouseWheelListener(wheelTracker);

	};

	public File getCurrentFile(){
		return this.currentFile;
	}

	public double getScale(){
		return this.scale;		
	}

	public void assignScale(double theScale){
		this.scale= theScale;
	}

	public void setCurrentFile(File currentFile){
		this.currentFile = currentFile;
	}

	public void setLookAndFeel(String clazz)
	{
		if (this != null)
		{
			try
			{
				UIManager.setLookAndFeel(clazz);
				SwingUtilities.updateComponentTreeUI(this);
			}	
			catch (Exception e1)
			{
				e1.printStackTrace();
			}
		}
	}

	public void setMenuBar(){
		//swing component
		menuBar = new JMenuBar();

		//filechooser
		fc = new JFileChooser();

		//build the file menu
		fileMenu = new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		fileMenu.getAccessibleContext().setAccessibleDescription(
				"The only menu in this program that has menu items");

		// A group of JMenuItems for exportMenu
		exprotMenu = new JMenu("Export to");
		menuItem1  = new JMenuItem("PNG file");
		menuItem1.addActionListener(new ExportPNGAction("export as","PNG",KeyEvent.VK_1));

		menuItem2  = new JMenuItem("GIF file");
		menuItem2.addActionListener(new ExportTIFAction("export as","GIF",KeyEvent.VK_2));

		exprotMenu.add(menuItem1);
		exprotMenu.add(menuItem2);

		fileMenu.add(exprotMenu);

		fileMenu.addSeparator();

		// Exit Menu for exitmneu 
		exitMenu = new JMenuItem("Exit");
		exitMenu.addActionListener(new ActionListener()
		{ public void actionPerformed(ActionEvent e)
		{
		   exit();
		}
		});			
		fileMenu.add(exitMenu);

		//build the view menu
		viewMenu = new JMenu("View");
		viewMenu.setMnemonic(KeyEvent.VK_V);
		viewMenu.getAccessibleContext().setAccessibleDescription(
				"The only menu in this program that has menu items");

		// A group of JMenuItems for layout
		layoutMenu = new JMenu("Layout");
                layoutMenu.setMnemonic(KeyEvent.VK_L);
		circleLayoutItem = new JMenuItem("Circle");
		hieraLayoutItem = new JMenuItem("Hierarchical");
		compactTreeLayoutItem = new JMenuItem("CompactTree");

		circleLayoutItem.addActionListener(new CircleLayoutAction());		
		hieraLayoutItem.addActionListener(new HieraLayoutAction());
		compactTreeLayoutItem.addActionListener(new CompactTreeLayoutAction());

		layoutMenu.add(circleLayoutItem);
		layoutMenu.add(hieraLayoutItem);
		layoutMenu.add(compactTreeLayoutItem);

		// A group of JMenuItems for ScaleMenu
		scaleMenu = new JMenu("Scale");
		menuItem1 = new JMenuItem("20%");
		menuItem2 = new JMenuItem("50%");
		menuItem3 = new JMenuItem("100%");

		menuItem1.addActionListener(new Scale20Action());
		menuItem2.addActionListener(new Scale50Action());
		menuItem3.addActionListener(new Scale100Action());

		scaleMenu.add(menuItem1);
		scaleMenu.add(menuItem2);
		scaleMenu.add(menuItem3);
		viewMenu.add(scaleMenu);

		// A group of JMenuItems for UndoMenu
		//UndoMenuItem = new JMenuItem("Undo",KeyEvent.VK_Z);
		//viewMenu.add(UndoMenuItem);

		//about menu
		helpMenu = new JMenu("Help");
		helpMenu.setMnemonic(KeyEvent.VK_H);
		aboutMenu = new JMenuItem("About GNV");
		aboutMenu.addActionListener(new ActionListener()
		{
			/*
			 * (non-Javadoc)
			 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
			 */
			public void actionPerformed(ActionEvent e)
			{
				about();
			}
		});

		helpMenu.add(aboutMenu);

		//	fileMenu.add(menuItem3);
		menuBar.add(fileMenu);
		menuBar.add(viewMenu);		
		menuBar.add(layoutMenu);
		menuBar.add(helpMenu);

		this.setJMenuBar(menuBar);
		this.setTitle(title);

		JFileChooser fc = new JFileChooser();	
	}

	public void createEdgeStyle(){ 
		Map<String, Object> stil = new HashMap<String, Object>();
		stil.put(mxConstants.STYLE_ROUNDED, true);
		stil.put(mxConstants.STYLE_EDGE, mxConstants.EDGESTYLE_ELBOW); 
		//	stil.put(mxConstants.STYLE_SHAPE, mxConstants.SHAPE_CONNECTOR);
		stil.put(mxConstants.STYLE_ENDARROW, mxConstants.ARROW_CLASSIC);
		stil.put(mxConstants.STYLE_ELBOW, mxConstants.ELBOW_VERTICAL);
		stil.put(mxConstants.STYLE_ALIGN, mxConstants.ALIGN_CENTER);
		stil.put(mxConstants.STYLE_STROKECOLOR, "#6482B9");
		stil.put(mxConstants.STYLE_FONTCOLOR, "#446299");
		stil.put(mxConstants.STYLE_ENDSIZE,10);
		stil.put(mxConstants.STYLE_STROKEWIDTH,2);
		this.graphStyleSheet.putCellStyle("RxnEdge",stil);
		//edgeStyle=elbowEdgeStyle;elbow=vertical;strokeWidth=3;strokeColor=#634E4E;endSize=3"
	}

	// compStyle = "shape=swimlane;fontSize=20;fontStyle=3;startSize=30;horizontal=false;autosize=1;align=top";
	public void createSwimLaneStyle(){ 
		Map<String, Object> stil = new HashMap<String, Object>();
		stil.put(mxConstants.STYLE_SHAPE, mxConstants.SHAPE_SWIMLANE);
		stil.put(mxConstants.STYLE_FONTSIZE, 40);
		stil.put(mxConstants.STYLE_FONTSTYLE,mxConstants.FONT_ITALIC);
		stil.put(mxConstants.STYLE_FONTCOLOR, "#446299");
		stil.put(mxConstants.STYLE_STARTSIZE,33);
		stil.put(mxConstants.STYLE_AUTOSIZE,1);
		stil.put(mxConstants.STYLE_HORIZONTAL, "false");
		this.graphStyleSheet.putCellStyle("SwimLane", stil);	
	}

	@SuppressWarnings("deprecation")
	public void setScale(double scale){
		// 
		this.assignScale(scale);	 

		//set edge style same as before scale
		this.graph.getView().setScale(scale);
		this.glycanCanvas.setGlycanScale(scale);
		this.createEdgeStyle();
		this.graph.setStylesheet(graphStyleSheet);
	}
        
        public void fit(){
                Dimension dimCustom = this.getSize();
            //    System.out.println("the size of frame is "+dimCustom.height+" "+dimCustom.width);
             //   System.out.println("the full size of the graph is "+defaultSize.height+" "+defaultSize.width);
                                
                double scaleY = ((double) dimCustom.height)/((double) defaultSize.height);
                double scaleX = ((double) dimCustom.width)/((double) defaultSize.width);
                
                double scaleToFit;
                double offset = 0.10;
                scaleToFit = Math.min(scaleX,scaleY)* (1-offset);
            //    System.out.println("scale factor is "+scaleToFit);
                if(scaleToFit<=0) scaleToFit=1;                
		this.setScale(scaleToFit);
        }

	public void fit(Dimension dimCustom){
		//System.out.println("dimension is "+defaultSize.getWidth()+defaultSize.getHeight());
                double scaleY = ((double) dimCustom.height)/((double) defaultSize.height);
                double scaleX = ((double) dimCustom.width)/((double) defaultSize.width);
                
                double scaleToFit=1;
                double offset = 0.05;
                scaleToFit = Math.min(scaleX,scaleY) - offset;
                System.out.println("scale factor is "+scaleToFit);
		this.setScale(scaleToFit);		
	}	

	private Dimension computeDefaultSize() {
		return null;
	}

	public GlycanNetGraphEditor(){
		this.swimlaneManager=new mxSwimlaneManager(graph);
		this.swimlaneManager.setAddEnabled(true);
		this.swimlaneManager.setResizeEnabled(true);
		this.swimlaneManager.setEnabled(true);
		
		this.theGlycanNets = new Vector<GlycanNet>();
		this.setBackground(Color.white);
		graphComponent = new mxGraphComponent(graph)
		{
			private static final long serialVersionUID = 4683716829748931448L;

			public GlycanSwingCanvas createCanvas()
			{
				return new GlycanSwingCanvas(this);							
			}
		};
		this.glycanCanvas = (GlycanSwingCanvas) graphComponent.getCanvas();
		this.theGraphicsOptions = this.glycanCanvas.getGraphicOptions();
                
		this.graphStyleSheet = this.graph.getStylesheet();

		mxIEventListener undoHandler = new mxIEventListener()
		{
			public void invoke(Object source, mxEventObject evt)
			{
				List<mxUndoableChange> changes = ((mxUndoableEdit) evt
						.getProperty("edit")).getChanges();
				graph.setSelectionCells(graph
						.getSelectionCellsForChanges(changes));
			}
		};

		// Updates the modified flag if the graph model changes
		this.graph.getModel().addListener(mxEvent.CHANGE, changeTracker);

		// Adds the command history to the model and view
		this.graph.getModel().addListener(mxEvent.UNDO, undoHandler);
		this.graph.getView().addListener(mxEvent.UNDO, undoHandler);

		// Keeps the selection in sync with the command history
		undoManager.addListener(mxEvent.UNDO, undoHandler);
		undoManager.addListener(mxEvent.REDO, undoHandler);		

		// Creates the graph outline component
		graphOutline = new mxGraphOutline(graphComponent);

		this.graphComponent.getGraphControl().addMouseListener(new MouseAdapter() {
			@Override
			public void mousePressed(MouseEvent e) {
			//	System.out.println("e.getButton() = " + e.getButton());
				if(e.getClickCount()==2){
                                Object cell = graphComponent.getCellAt(e.getX(), e.getY());
				//System.out.println("Mouse click in graph component");
				if (cell != null) {
					try{
						GlycanSpecies theGlycanSpecies = (GlycanSpecies) ((mxCell) ((mxCell) cell).getValue()).getValue();
					//	System.out.println("the glycan name is "+theGlycanSpecies.getName());
						try {
						    String dbName = "glycomedb";	
                                                    linkToGlycomeDBsites(theGlycanSpecies,dbName);
						} catch (IOException ex) {
							System.out.println(ex.getMessage());
							//System.out.println();
						}
					}catch(Exception ee) {
						//
					}
					//System.out.println("cell=" + graph.getLabel(cell));
				}
                                }  
			}

            public void linkToGlycomeDBsites(GlycanSpecies theGlycanSpecies,String dbName) throws URISyntaxException, IOException {
                //   Process pc = Runtime.getRuntime().exec("cmd.exe /c start http://www.genome.jp/dbget-bin/www_bget?gl:G00109");

                // String keggIDFromCell = theGlycanSpecies.getDbIDs().get("KEGG");
                String glycomeID = theGlycanSpecies.getDbIDs().get(dbName); 
                String glycanURI;
                if((dbName=="glycomedb")&&(!glycomeID.equalsIgnoreCase("-1")))
                    glycanURI = prefix_GLYCOMEDB + glycomeID;
                else if((dbName=="KEGG")&&(!glycomeID.equalsIgnoreCase("-1")))
                    glycanURI = prefix_KEGG + glycomeID;
                else 
                    glycanURI = link_GLYCOMEDB;
                
                URI uri = new java.net.URI(glycanURI);
                java.awt.Desktop desktop = java.awt.Desktop.getDesktop();
                desktop.browse(uri);
            }
		});		
	}

	public GlycanNetGraphEditor(String sbmlFileName) throws Exception{
		//	String SBMLfileName="UB1997Model.xml";
		this.swimlaneManager=new mxSwimlaneManager(this.graph);
		this.swimlaneManager.setAddEnabled(true);
		this.swimlaneManager.setResizeEnabled(true);
		this.swimlaneManager.setEnabled(true);
		GlycanNetRead ubNetModel = new GlycanNetRead(sbmlFileName);
		//ubNetModel.setVisible(false);
		this.theGlycanNets = ubNetModel.theGlycanRxnNets;

		this.setBackground(Color.white);
		graphComponent = new mxGraphComponent(graph)
		{
			private static final long serialVersionUID = 4683716829748931448L;

			public GlycanSwingCanvas createCanvas()
			{
			   return new GlycanSwingCanvas(this);
			}
		};
		this.glycanCanvas = (GlycanSwingCanvas) graphComponent.getCanvas();
		this.theGraphicsOptions = this.glycanCanvas.getGraphicOptions();
		this.graphStyleSheet = this.graph.getStylesheet();                                                                 
	}

	public void visNet(Vector<GlycanNet> theNetworks) throws Exception{
		double left=10;
		double top=10;
		double width = 400;
		double modelheight = 400;

		graph.getModel().beginUpdate();
		if(theNetworks.isEmpty()) return;	    
		try {
			Object parent = graph.getDefaultParent();
			double nCompts = theNetworks.size();
			double factionInt = 0.1;
			double height =  modelheight/(nCompts+(nCompts-1)*factionInt);
			double interval = factionInt*height;

			for(int i=0; i<theNetworks.size();i++) {
				//System.out.println(" "+i+" 's compartmtent");
				Object comptCell = this.graph.insertVertex(parent,null, 
						theNetworks.get(i).getCompartName(),left,top,width,height,
						"SwimLane");
				this.glycanNetVertexManager.put(theNetworks.get(i),comptCell);
				top = top + height + interval;
				visNet(theNetworks.get(i),comptCell);
				//System.out.println("manager "+this.swimlaneManager.isAddEnabled());	
			}		
		}finally{
			graph.getModel().endUpdate();
		}
	}

	public void visNet(GlycanNet theNetwork) throws Exception {
		double left=0;
		double top=0;
		double width =  200;
		double height = 200;
		//System.out.println("the compartment name is "+theNetwork.getCompartName());

		graph.getModel().beginUpdate();	
		try{
			mxCell comptCell = (mxCell) this.graph.insertVertex(graph.getDefaultParent(),null,
					theNetwork.getCompartName(),left,top,width,height,//,"SwimLane");
			"shape=swimlane;fontSize=20;fontStyle=3;startSize=30;horizontal=false;autosize=1;");
			//comptCell.setId(theNetwork.getCompartName());
			this.glycanNetVertexManager.put(theNetwork,comptCell);	
			visNet(theNetwork,comptCell);	
		}
		finally
		{
			graph.getModel().endUpdate();
		}
		graph.repaint();

	}

	public void visNet(GlycanNet theNetwork,Object compt) throws Exception {

		//species
		Vector<GlycanSpecies> theSpecies = theNetwork.getStructures();
		Vector<GlycanRxn> theRxns = theNetwork.getRxns();	 		

		//add species to the graph			
		addGlycansToGraph(compt,new mxPoint(10,10),theSpecies);

		//add reactions to the graph
		addRxnsToGraph(compt,theRxns);
		compt = this.graph.updateCellSize(compt,false);	
	}

	private void addRxnsToGraph(Object compt, Vector<GlycanRxn> theRxns) throws Exception {
		for(int i=0;i<theRxns.size();i++) addRxnToGraph(compt,theRxns.get(i));
	}

	private void addRxnToGraph(Object compt,
			GlycanRxn glycanRxn) throws Exception {
		Object reac =  glycanVertexManager.get(glycanRxn.getReactants().get(0));
		Object prod =  glycanVertexManager.get(glycanRxn.getProducts().get(0));
		if((reac==null)||(prod==null)) {
                    throw new Exception("Current versiono only supports 1->1 reaction");
                    //System.out.println("rxn is not added successfully");
                }
		mxCell glycanEdge = (mxCell) graph.insertEdge(compt, null, glycanRxn.getEnz(), reac, prod, "RxnEdge");
		glycanRxnEdgeManager.put(glycanRxn, glycanEdge);		
	}


	private void addGlycansToGraph(Object compt, mxPoint left_top,
			Vector<GlycanSpecies> theSpecies) {
		for(int i=0;i<theSpecies.size();i++)  
			left_top = addGlycanToGraph(compt,left_top,theSpecies.get(i));
	}

	private mxPoint addGlycanToGraph(Object compt, mxPoint left_top,
			GlycanSpecies glycanSpecies) {		
		GlycanSpecies glycanToVis =   glycanSpecies;
		int left_posX =  (int) left_top.getX();
		int left_posY =  (int) left_top.getY()+30;

		Rectangle glycanBox = this.glycanCanvas.getGlycanRenderer().computeBoundingBoxes(glycanToVis.getStructure(), 
				left_posX, left_posY, true, true,this.glycanCanvas.thePosManager,this.glycanCanvas.theBBoxManager);
		//	System.out.println("the size of the box is "+glycanBox.getX()+" "+glycanBox.getY()+" width: "+glycanBox.getWidth()+" height: "+glycanBox.getHeight());

		mxCell cellGlycan = new mxCell(glycanToVis);		
		Object vertexGlycan = graph.insertVertex(compt, null, cellGlycan, glycanBox.getX(), glycanBox.getY(), glycanBox.getWidth(),glycanBox.getHeight());
		glycanVertexManager.put(glycanToVis,vertexGlycan);		
		return new mxPoint(glycanBox.getMinX(),glycanBox.getMaxY());
	}

	@SuppressWarnings("deprecation")
	public String setSwimLaneStype(){
		String rootStyle = "strokeColor=black";
		rootStyle = mxStyleUtils.setStyle(rootStyle, mxConstants.STYLE_SHAPE, mxConstants.SHAPE_SWIMLANE);
		//,mxConstants.FONT_BOLD,mxConstants.DEFAULT_FONTSIZE);
		rootStyle = mxStyleUtils.setStyleFlag(rootStyle, mxConstants.STYLE_FONTSIZE,20, true);
		rootStyle = mxStyleUtils.setStyleFlag(rootStyle, mxConstants.STYLE_FONTSTYLE,mxConstants.FONT_BOLD, true);
		rootStyle = mxStyleUtils.setStyleFlag(rootStyle, mxConstants.STYLE_ENDSIZE,30, true);
		rootStyle = mxStyleUtils.setStyleFlag(rootStyle, mxConstants.STYLE_AUTOSIZE,1, true);
		rootStyle = mxStyleUtils.setStyle(rootStyle, mxConstants.ELBOW_VERTICAL,"true");
		return rootStyle;
	}

	/*
	 * 
	 * 
	 */
	public void visNetFromSBMLFile(String SBMLfileName,String glycanFileFormat,String frameDimensionChoice,
        Dimension frameDimension,boolean showRedEnd, boolean showMass, boolean showLinkage,String netLayout) throws Exception
			{
                this.visNetFromSBMLFile(SBMLfileName,glycanFileFormat,frameDimensionChoice,frameDimension,showRedEnd,showMass,showLinkage,netLayout,null,false);
}
        
        public void visNetFromSBMLFile(String SBMLfileName,String glycanFileFormat,String frameDimensionChoice,
        Dimension frameDimension,boolean showRedEnd, boolean showMass, boolean showLinkage,String netLayout,HashMap<String,Integer> glycomedbIDs) throws Exception
			{
                this.visNetFromSBMLFile(SBMLfileName,glycanFileFormat,frameDimensionChoice,frameDimension,showRedEnd,showMass,showLinkage,netLayout,glycomedbIDs,false);
        }
        
	public void visNetFromSBMLFile(String SBMLfileName,String glycanFileFormat,String frameDimensionChoice, Dimension frameDimension, 
			boolean showRedEnd, boolean showMass, boolean showLinkage,String netLayout,HashMap<String,Integer> glycomedbIDs,boolean isfit) throws Exception
			{
                this.setFrameOptions(frameDimensionChoice, frameDimension, showRedEnd, showMass, showLinkage);                            
                if(glycanFileFormat==null) glycanFileFormat="glycoct_xml";
		this.theGlycanNets = this.readSBMLFile(SBMLfileName,glycanFileFormat);
                if(glycomedbIDs!=null){
                    this.assignGlycomeID(this.theGlycanNets,glycomedbIDs);
                }
            //    System.out.println("the size of the network is: debugging point 1"+this.theGlycanNets.size());
                
		//this.visNet(firstCompt);
		this.visNet(this.theGlycanNets);

		//frame.visNet(frame.theGlycanNets);
		new mxRubberband(this.graphComponent);
                this.layoutGlycanNetGraph(netLayout,isfit);
                this.setVisible(true);
		
	//	System.out.println("size of the network is "+glycanNetVertexManager.size());
	//	System.out.println("the default size for the network is " + 
         //       this.defaultSize.getHeight() + " " + this.defaultSize.getWidth());
	//	System.out.println("done for visulization");		
}

	private void assignGlycomeID(Vector<GlycanNet> theGlycanNets,
			HashMap<String,Integer> glycomedbIDs) {
		for(int i=0;i<theGlycanNets.size();i++){
			GlycanNet theithGlycanNet = theGlycanNets.get(i);
			assignGlycomeID(theithGlycanNet,glycomedbIDs);
		}		
	}

	private void assignGlycomeID(GlycanNet theithGlycanNet, HashMap<String,Integer> glycomedbIDs) {
		// TODO Auto-generated method stub
		for(int i=0;i<theithGlycanNet.getStructures().size();i++){
			GlycanSpecies theGlycan = theithGlycanNet.getStructures().get(i);
			String glycoctGlycanStruct = theGlycan.getStructure().toGlycoCT();
			theGlycan.addDbID("glycomedb", Integer.toString(glycomedbIDs.get(glycoctGlycanStruct)));
		}
	}

	public void layoutGlycanNetGraph(String netLayout){
                this.layoutGlycanNetGraph(netLayout,false);
	}       

	private void setDefaultSize(int width,int height) {
		// TODO Auto-generated method stub
		this.defaultSize= new Dimension(width,height);
              //  System.out.println("the default size is "+width+" "+height);
	}

	private Dimension computePreferredSize() {
		// TODO Auto-generated method stub
		return null;
	}

	public Vector<GlycanNet> readSBMLFile(String SBMLfileName) throws Exception {
		String defaultGlycanFileFormat ="glycoct_xml";
		return readSBMLFile(SBMLfileName,defaultGlycanFileFormat);		
	}

	public Vector<GlycanNet> readSBMLFile(String SBMLFileName,String glycanFileFormat) throws Exception{
		//String SBMLfileName="UB1997Model.xml";
		glycanNetModel = new GlycanNetRead();
		glycanNetModel.readSBMLfromFile(SBMLFileName);		
                return glycanNetModel.theGlycanRxnNets;
	}

	protected mxIGraphLayout createMxGraphLayout(String layoutformatString,mxGraph graph){
		if(layoutformatString.equalsIgnoreCase("Hierarchical")){
			mxHierarchicalLayout layout = new mxHierarchicalLayout(graph,SwingConstants.NORTH); 
			return layout;
		} else if(layoutformatString.equalsIgnoreCase("Orthogonal"))
			return new mxOrthogonalLayout(graph);
		else if(layoutformatString.equalsIgnoreCase("Circle")){
			mxCircleLayout layoutCircle = new mxCircleLayout(graph);
			//layoutCircle.setResetEdges(true);
			//layoutCircle.setDisableEdgeStyle(true);
			return layoutCircle;
		}else if(layoutformatString.equalsIgnoreCase("CompactTree"))
			return new mxCompactTreeLayout(graph);
		else if(layoutformatString.equalsIgnoreCase("EdgeLable"))
			return new mxEdgeLabelLayout(graph);
		else if(layoutformatString.equalsIgnoreCase("Organic"))
			return new mxOrganicLayout(graph);
		else if(layoutformatString.equalsIgnoreCase("ParallelEdge"))
			return new mxParallelEdgeLayout(graph);
		else if(layoutformatString.equalsIgnoreCase("Partition"))
			return new mxPartitionLayout(graph);
		else if(layoutformatString.equalsIgnoreCase("Stack"))
			return new mxStackLayout(graph);
		else
			return new mxCompactTreeLayout(graph);
	}

	protected void layoutGraph(Object objectToLayout, mxIGraphLayout themxIGraphLayout){
		themxIGraphLayout.execute(objectToLayout);
	}

	/**
	 * 
	 */
	protected void mouseWheelMoved(MouseWheelEvent e)
	{
		if (e.getWheelRotation() < 0)
		{
			this.graphComponent.zoomIn();
		}
		else
		{
			this.graphComponent.zoomOut();
		}

		//zoom in glycan structures
		this.setScale(graphComponent.getGraph().getView().getScale());
		/*status(mxResources.get("scale") + ": "
					+ (int) (100 * graphComponent.getGraph().getView().getScale())
					+ "%");*/
	}

	protected void export(String filename){
		try
		{
			String ext = filename
					.substring(filename.lastIndexOf('.') + 1);

			if (ext.equalsIgnoreCase("svg"))
			{
				mxSvgCanvas canvas = (mxSvgCanvas) mxCellRenderer
						.drawCells(this.graph, null, 1, null,
								new CanvasFactory()
						{
							public mxICanvas createCanvas(
									int width, int height)
							{
								mxSvgCanvas canvas = new mxSvgCanvas(
										mxDomUtils.createSvgDocument(
												width, height));
								canvas.setEmbedded(true);

								return canvas;
							}

						});

				mxUtils.writeFile(mxXmlUtils.getXml(canvas.getDocument()),
						filename);
			}
			else if (ext.equalsIgnoreCase("html"))
			{
				mxUtils.writeFile(mxXmlUtils.getXml(mxCellRenderer
						.createHtmlDocument(graph, null, 1, null, null)
						.getDocumentElement()), filename);
			}
			else
			{
				Color bg = null;

				if ((!ext.equalsIgnoreCase("gif") && !ext
						.equalsIgnoreCase("png")))
				{
					bg = graphComponent.getBackground();
				}
				else
				{
					//double scaleToDraw = 1;
					//this.setScale(scaleToDraw);
					BufferedImage image = mxCellRenderer
							.createBufferedImage(graph, null, this.getScale(), Color.WHITE,
									graphComponent.isAntiAlias(), null,graphComponent.getCanvas());
				//	System.out.println("debugging in export ");

					if (image != null)
					{
						ImageIO.write(image, ext, new File(filename));
					}					
				}
			}
		}
		catch (Throwable ex)
		{
			ex.printStackTrace();
			JOptionPane.showMessageDialog(graphComponent,
					ex.toString(), mxResources.get("error"),
					JOptionPane.ERROR_MESSAGE);
		}
	}

    public void layoutGlycanNetGraph(String netLayout, boolean isfit) {
                this.graph.setAutoOrigin(true);
		this.graph.setAutoSizeCells(true);	
		this.graph.setCellsEditable(true);
		this.graph.setCellsResizable(true);
            
		double top=0;
                mxCell comptObj = (mxCell) glycanNetVertexManager.get(this.theGlycanNets.get(0));			
		this.layoutGraph(comptObj,this.createMxGraphLayout(netLayout,this.graph));
		comptObj.getGeometry().setY(top);
		comptObj = (mxCell) this.graph.updateCellSize(comptObj, false);
		double height = comptObj.getGeometry().getHeight();
                double width = comptObj.getGeometry().getWidth();
           //     System.out.println("width is "+width);
                int defaultSizeX = (int) comptObj.getGeometry().getWidth();
		//System.out.println("top of the cell is "+comptObj.getGeometry().getHeight()+" new height is  "+height);
                // System.out.println("width is "+defaultSizeX);
		for(int i=1;i<theGlycanNets.size();i++){			
			comptObj = (mxCell) glycanNetVertexManager.get(theGlycanNets.get(i));			
			this.layoutGraph(comptObj,this.createMxGraphLayout(netLayout,this.graph));
                        top = top + height;
			comptObj.getGeometry().setY(top);
                        comptObj.getGeometry().setHeight(height);
			//comptObj = (mxCell) this.graph.updateCellSize(comptObj, false);     		        
		//	System.out.println("top of the cell is "+comptObj.getGeometry().getHeight()+" new height is  "+height);
		}
                
                int defaultSizeY = (int) (top+height);
                this.setDefaultSize(defaultSizeX,defaultSizeY);
              //  System.out.println("the default size in layout is "+defaultSizeX+" "+defaultSizeY);
		this.graph.repaint();
                if(isfit) fit();
    
    }

	public class ExportPNGAction extends AbstractAction{
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public ExportPNGAction(String text,String desc,
				Integer mnemonic){
			super(text);
			putValue("Short_Descriptioin",desc);
			putValue("Mnemonic Key",mnemonic);
		}

		@Override
		public void actionPerformed(ActionEvent e) {
			if(fc == null){
				fc = new JFileChooser();
			}

			//Add a custom file filter and disable the default file filter
			FileFilter ft = new FileNameExtensionFilter("PNG File","PNG");
			fc.addChoosableFileFilter(ft);
			fc.setAcceptAllFileFilterUsed(false);

			int returnVal = fc.showSaveDialog(GlycanNetGraphEditor.this);
			if(returnVal == JFileChooser.APPROVE_OPTION) {
				File chosenFile = fc.getSelectedFile();
				try {
					System.out.println(chosenFile.getCanonicalPath());
					String newFileName =chosenFile.getCanonicalPath();
					if(!newFileName.endsWith("png")){
						newFileName=newFileName+".png";
					}
					System.out.println(newFileName);
					export(newFileName);
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}			
		}		
	}

	public class Scale20Action extends AbstractAction{

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			setScale(0.2);	
		}	
	}

	public class Scale50Action extends AbstractAction{

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			setScale(0.5);	
		}	
	}

	public class Scale100Action extends AbstractAction{

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			setScale(1);	
		}	
	}



	public class ExportTIFAction extends AbstractAction{
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public ExportTIFAction(String text,String desc,
				Integer mnemonic){
			super(text);
			putValue("Short_Descriptioin",desc);
			putValue("Mnemonic Key",mnemonic);
		}

		@Override
		public void actionPerformed(ActionEvent e) {
			if(fc == null){
				fc = new JFileChooser();
			}

			//Add a custom file filter and disable the default file filter
			FileFilter ft = new FileNameExtensionFilter("gif File","gif");
			fc.addChoosableFileFilter(ft);
			fc.setAcceptAllFileFilterUsed(false);

			int returnVal = fc.showSaveDialog(GlycanNetGraphEditor.this);
			if(returnVal == JFileChooser.APPROVE_OPTION) {
				File chosenFile = fc.getSelectedFile();
				try {
					System.out.println(chosenFile.getCanonicalPath());
					String newFileName =chosenFile.getCanonicalPath();
					if(!newFileName.endsWith("gif")){
						newFileName=newFileName+".gif";
					}
					System.out.println(newFileName);
					export(newFileName);
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}			
		}		
	}

	public void exit(){
		// TODO Auto-generated method stub
		if(this!=null){
			this.dispose();
		}		
	}

	public mxGraphComponent getGraphComponent() {
		// TODO Auto-generated method stub
		return this.graphComponent;
	}

	public Object getUndoManager() {
		// TODO Auto-generated method stub
		return null;
	}

	public class CircleLayoutAction extends AbstractAction{

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		@Override
		public void actionPerformed(ActionEvent e) {
		//	System.out.println("the Circle");
			layoutGlycanNetGraph("Circle");

			//Dimension frameDimn = new Dimension(getWidth(),getHeight());
			//fit(frameDimn);
		}	
	}

	public class HieraLayoutAction extends AbstractAction{

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		@Override
		public void actionPerformed(ActionEvent e) {
			//System.out.println("the hier");
			layoutGlycanNetGraph("Hierarchical"); 

			//layoutManager.destroy();
			//setLayoutManager("Hierarchical");
			//Dimension frameDimn = new Dimension(getWidth(),getHeight());
			//fit(frameDimn);
		}	
	}

	public class CompactTreeLayoutAction extends AbstractAction{

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		@Override
		public void actionPerformed(ActionEvent e){
			//System.out.println("the compact tree");
			layoutGlycanNetGraph("CompactTree");

			//layoutManager.destroy();
			//setLayoutManager("CompactTree");
			//Dimension frameDimn = new Dimension(getWidth(),getHeight());
			//fit(frameDimn);
		}	
	}

	public void about()
	{
		if (this != null)
		{
			AboutFrame about = new AboutFrame(this);
			about.setModal(true);

			// Centers inside the application frame
			int x = this.getX() + (this.getWidth() - about.getWidth()) / 2;
			int y = this.getY() + (this.getHeight() - about.getHeight()) / 2;
			about.setLocation(x, y);
			//System.out.println("call about");

			// Shows the modal dialog and waits
			about.setVisible(true);
		}
	}
}



