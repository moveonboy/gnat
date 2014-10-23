package org.glyco;
/**
 * 
 */

/**
 * @author gangliu
 *
 */
import java.awt.BorderLayout;
import java.net.MalformedURLException;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import org.eurocarbdb.application.glycanbuilder.*;

public class GlycanRead extends JFrame{

	/**
	 * 
	 */
    private static final long serialVersionUID = -204730955259211308L;
    public BuilderWorkspace theWorkspace;
    public GlycanDocument theDoc;
    public GlycanCanvas theCanvas;
    public JScrollPane sp;

	/**
	 * 
	 */
	//private static final long serialVersionUID = 1L;
    public GlycanRead(){
      	            theWorkspace = new BuilderWorkspace();
		    theWorkspace.setAutoSave(true);
		    GraphicOptions theGraphicOptions = theWorkspace.getGraphicOptions();
		    theGraphicOptions.SHOW_MASSES_CANVAS=false;
		    theGraphicOptions.SHOW_REDEND_CANVAS=false;
		    theGraphicOptions.SHOW_INFO=true;
		    theWorkspace.setDisplay("compact");
		    
		    // create singletons 
		    this.theDoc = theWorkspace.getStructures();
		    
		    // set the layout
		    this.getContentPane().setLayout(new BorderLayout());	   
						    	    
			// create canvas
		    this.theCanvas = new GlycanCanvas(this,theWorkspace);
			 
			 // set the canvas
			 this.sp = new JScrollPane(theCanvas);
			 theCanvas.setScrollPane(sp);
			 this.getContentPane().add(sp,BorderLayout.CENTER);
			 //repaint();
			 
			// setIconImage(FileUtils.defaultThemeManager.getImageIcon("logo").getImage());
		     setSize(400, 300);        
		     setLocationRelativeTo(null);
    }
    
    public void read(Glycan theGlycan) {
    	theDoc.addStructure(theGlycan);
    }
	
	public void read(String fileName,String format) throws MalformedURLException {
	    //read structure 	   	
		theDoc.importFrom(fileName, format);				  
	}
	
	public void readFromString(String glycanString,String format){
		theDoc.importFromString(glycanString, format);
	}
	
	public void readFromString(String glycanStructStr) throws MalformedURLException {
	 
		 //read structure 	
		theDoc.importFromString(glycanStructStr, "glycoct_xml");
    }
	
	public void read() throws MalformedURLException {
	         // create the default workspace
	         String fileName="highmannose_linucs.xml";
		 String format = "linucs";
		 this.theDoc.importFrom(fileName, format);
   }
	
	public void exportToXMLFile(String outputFileName,String format){
		if(!this.theDoc.isEmpty()){
		  this.theDoc.exportTo(outputFileName, format);
		} else {
		  System.out.println("Please provide glycan structures");
		}
	}
	
	public String exportToString(String format){
		if(!this.theDoc.isEmpty()){
		  return this.theDoc.toString(format);
		} else {
		  System.out.println("Please provide glycan structures");
		}
		return null;
	}
	

	/**
	 * @param args
	 * @throws MalformedURLException 
	 */
	public static void main(String[] args) throws MalformedURLException {
		//test usage of glycan builder
	    String fileName="highmannose_linucs.xml";
	    System.out.println("testing new output format");
        // String format = "linucs";
        String format = "iupac_condenced";
        String linucsformat = "linucs";
        String cfgformat = "cfg";
        String glycanIUPACString = " Glc a#Sp8";
        String glycanLinucsString = "[][b-D-GlcpNAc]{[(4+1)][b-D-GlcpNAc]{[(4+1)][b-D-Manp]{[(3+1)][a-D-Manp]{[(2+1)][a-D-Manp]{}}[(6+1)][a-D-Manp]{[(3+1)][a-D-Manp]{}[(6+1)][a-D-Manp]{}}}}}";
        String glycanCFGString = "Ub3Ab3Ab4Xb3;S";
		GlycanRead defaultGlycan= new GlycanRead();
		defaultGlycan.setVisible(true);
		//defaultGlycan.read(fileName,format);
		defaultGlycan.readFromString(glycanIUPACString,cfgformat);
		String iupcaString = defaultGlycan.exportToString("iupac_condenced");
		System.out.println("new glyan string is "+iupcaString);
		
		
	//	defaultGlycan.exportToXMLFile("highmannose_linucs3.xml", "linucs");
	    System.out.println("succeed");
	}

}
