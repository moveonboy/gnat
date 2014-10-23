package org.glyco;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Vector;
import javax.xml.stream.XMLStreamException;
import org.eurocarbdb.application.glycanbuilder.BuilderWorkspace;
import org.eurocarbdb.application.glycanbuilder.Glycan;
import org.eurocarbdb.application.glycanbuilder.GlycanDocument;
import org.sbml.jsbml.*;


public class GlycanNetRead{
	/*
	 *  
	 */
	private static final long serialVersionUID = 1L;
	public Vector<GlycanNet> theGlycanRxnNets;
        public Vector<GlycanSpecies> theGlycanSpecies;
        
        public HashMap<String,GlycanSpecies> glycanIndexMap = new HashMap<String,GlycanSpecies>();
        public HashMap<GlycanNet,String> glycanNetsIndexMap = new HashMap<GlycanNet,String>();
        
        public Vector<GlycanRxn> theGlycanRxns;
        public Vector<Glycan> theGlycanArray; 
        public Vector<String> theCompartments;        
        
	public GlycanNetRead(){
		this.theGlycanRxnNets = new Vector<GlycanNet>();                
                this.theGlycanRxns    = new Vector<GlycanRxn>();
                this.theGlycanSpecies = new Vector<GlycanSpecies>();
                this.theGlycanArray   = new Vector<Glycan>();
                this.theCompartments  = new Vector<String>();
	}
        
        /*
         * read the GlycanNets from SBML
         */
	public GlycanNetRead(String sbmlFileName) throws Exception {
		this();
		readSBMLfromFile(sbmlFileName);
        }
        
        public GlycanNetRead(Vector<GlycanNet> theGlycanRxnNets)throws Exception {
		this();
                this.theGlycanRxnNets = theGlycanRxnNets;
                this.setGlycanSpeciesRxnCompts(theGlycanRxnNets);                                
        }
        
        public void setGlycanSpeciesRxnCompts(Vector<GlycanNet> theGlycanNets){
            for(int i=0;i<theGlycanNets.size();i++){
               GlycanNet theNet = theGlycanNets.get(i);
               //add compartment name
               theCompartments.add(theNet.getCompartName())  ;
               
               //add glycan
               for(int j=0;j<theNet.getStructures().size();j++){
                  theGlycanArray.add(theNet.getStructures().get(j).getStructure());
                  theGlycanSpecies.add(theNet.getStructures().get(j));
                  glycanIndexMap.put(theNet.getStructures().get(j).getStructure().toGlycoCT(),
                          theNet.getStructures().get(j));                  
               }
               
               //add glycan rxn
               for(int j=0;j<theNet.getRxns().size();j++){
                  this.theGlycanRxns.add(theNet.getRxns().get(j));                  
               }               
            }
        }
        
   //     static{
   //		System.loadLibrary("sbmlj");
   //	}

	/**
	 * @param args
	 * @throws Exception 
	 * @throws XMLStreamException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws Exception 
	{
		//test usage of glycan builder
		String SBMLfileName="UB1997Model.xml";
		GlycanNetRead ubNetModel = new GlycanNetRead();
		
		// read sbml file 
		System.out.println("before sbml import, the size of the network is "+ubNetModel.theGlycanRxnNets.size());
		
		ubNetModel.readSBMLfromFile(SBMLfileName);//
                ubNetModel.theGlycanRxnNets.remove(3);
//		ubNetModel.theGlyNet.remove(2);
//		ubNetModel.theGlyNet.remove(1);		
		System.out.println("after sbml import, the size of the network is "+ubNetModel.theGlycanRxnNets.size());		
		System.out.println("succeed");
	}
        
        public void readStructures(Model theSBMLModel) throws MalformedURLException {
		// this.glycanArray.setVisible(true);
		GlycanDocument theGlycanDocument = new GlycanDocument(new BuilderWorkspace());
				
		for(int i = 0; i < (int) theSBMLModel.getNumSpecies() ; i++)
		{
			 Species sp = theSBMLModel.getSpecies(i);
			// String a = sp.getNotesString();
			// a = note2GlycanStr(a);
                         String a=sp.getAnnotationString();
                         a = annotation2GlycanStr(a);
			 theGlycanDocument.importFromString(a, "glycoct_xml");
		}
		
		for(int i = 0; i < (int) theSBMLModel.getNumSpecies(); i++)
		{
		    Glycan theGlycan = theGlycanDocument.getStructures().get(i);
                    String speciesName;
                    speciesName = theSBMLModel.getSpecies(i).getName();                    
                    GlycanSpecies newSpecies = new GlycanSpecies(theGlycan,
                            speciesName,theSBMLModel.getSpecies(i).getCompartment());
                    newSpecies.setID(theSBMLModel.getSpecies(i).getId());
                    this.theGlycanSpecies.add(newSpecies);
                    this.theGlycanArray.add(theGlycan);
		    this.glycanIndexMap.put(theGlycan.toGlycoCT(),newSpecies);			 
		}
                
         /*       for(GlycanSpecies oneGlycanSpecies: this.theGlycanSpecies){
                   System.out.println("species "+" = "+oneGlycanSpecies.getName());
                }                
          */
	}
        
        public void readRxns(Model theSBMLModel) throws Exception {
             for(int i=0;i<(int) theSBMLModel.getNumReactions();i++){
		    GlycanRxn theGlycanRxn = sbml2GlycanRxn(theSBMLModel.getReaction(i));
                    this.theGlycanRxns.add(theGlycanRxn);
             }
        }
        
       public GlycanRxn sbml2GlycanRxn(Reaction sbmlRxn) 
               throws Exception {
		
			GlycanRxn theGlycanRxn;
                        GlycanSpecies reac=null,prod=null;
			
                        ListOf<SpeciesReference> reactRefList = sbmlRxn.getListOfReactants();
                        ListOf<SpeciesReference> prodRefList =  sbmlRxn.getListOfProducts();
                        
                        if(reactRefList.size()> 1 || prodRefList.size()>1){
                           throw new Exception("Current reaction format "
                                   + "only supprt one reaction and one product"); 
                        }
                        
                        if(reactRefList.size()== 1 && reactRefList.get(0).getSpecies()!= null){ 
                            for(GlycanSpecies oneGlycanSpecies: this.theGlycanSpecies){
                                if( (reactRefList.get(0).getSpecies().equalsIgnoreCase(oneGlycanSpecies.getName())) ||
                                    (reactRefList.get(0).getSpecies().equalsIgnoreCase(oneGlycanSpecies.getID())))   
                                       {
                                    reac = oneGlycanSpecies;
                                }
                            }
                        } else {
                           reac = null;
                        } 
                        
                        if(prodRefList.size()== 1 && prodRefList.get(0).getSpecies()!= null){
                            for(GlycanSpecies oneGlycanSpecies: this.theGlycanSpecies){
                                if((oneGlycanSpecies.getName().equalsIgnoreCase(prodRefList.get(0).getSpecies())) ||
                                     (prodRefList.get(0).getSpecies().equalsIgnoreCase(oneGlycanSpecies.getID())))  
                                    prod = oneGlycanSpecies;                                
                            }
                        } else {
                            prod = null;
                        }         
                        			
                       	if((reac==null)&&(prod==null)) 
                        {
                     //       System.out.println( this.theGlycanSpecies.size());
                     /*       for(GlycanSpecies oneGlycanSpecies: this.theGlycanSpecies){
                                 System.out.println(oneGlycanSpecies.getName());
                            }
                            * 
                            */
                            
                     /*       for(int i=0;i<this.theGlycanSpecies.size();i++){
                                  System.out.println("species "+i+" = "+this.theGlycanSpecies.get(i).getName());
                            }
                            * 
                            */
                            
                 /*           System.out.println("reaction is "+reactRefList.get(0).getSpecies()
                                    +" "+prodRefList.get(0).getSpecies());
                                    * 
                                    */
                            throw new Exception("reaction set up is not right");
                        }
                        
                        theGlycanRxn = new GlycanRxn(reac,prod);
			
			return theGlycanRxn;
	}
       
        public void readCompartments(Model theSBMLModel){
            for(int i=0;i<theSBMLModel.getNumCompartments();i++){
                String comptName;
                if(theSBMLModel.getCompartment(i).getName().isEmpty())
                    comptName = theSBMLModel.getCompartment(i).getId();
                else
                    comptName = theSBMLModel.getCompartment(i).getName();
                
                this.theCompartments.add(comptName);
            }
        }
        
        public void readPathwayInCompts(Model ubModel) throws MalformedURLException, Exception{
              this.readPathwayInCompts(ubModel,true);
        }
        
        public void readPathwayInCompts() throws MalformedURLException, Exception{
              this.readPathwayInCompts(true);   
        }
        
        public void readPathwayInCompts(boolean reset) throws MalformedURLException, Exception{
                int numCompartments =(int) this.theCompartments.size();
		if(reset){
			this.theGlycanRxnNets.clear();
		}			
		
                //create the compartments and add the compartment names
		Vector<GlycanNet> theCompts = new Vector();
                for(int j=0;j<numCompartments;j++){
		    GlycanNet theGlycanCompt = new GlycanNet();
		    theGlycanCompt.setCompartName(this.theCompartments.get(j));
                    theCompts.add(theGlycanCompt);
               }
                
		//assign the species to each compartment
                for(int i = 0; i < (int) this.theGlycanSpecies.size() ; i++)
		{
		    for(GlycanNet theGlycanCompt:theCompts){
                      if(theGlycanSpecies.get(i).getCompartname().equalsIgnoreCase(theGlycanCompt.getCompartName()))
                      {
			theGlycanCompt.addStructure(theGlycanSpecies.get(i));
                        break;
                      }
                    }                    
		}
                
                //add the reactions to each compartment based on compartment information
               for(GlycanNet theGlycanCompt:theCompts){
           //          System.out.println("name of the compartment is "+theGlycanCompt.getCompartName());		
           //          System.out.println("number of the species is "+theGlycanCompt.getStructures().size());
			
		     //add rxns to each single compartment
                        int numInFlux=0,numOutFlux=0,numAdded=0,numNotAdded = 0;
                        
			for(int i=0;i<this.theGlycanRxns.size();i++){
                            boolean isProductsinOneCompt  = true;
                            boolean isReactantsinOneCompt = true; 
			    
                            GlycanRxn glycanRxn_compt = this.theGlycanRxns.get(i);
                                                            
                            if( (glycanRxn_compt.getType().equalsIgnoreCase("reaction")))
                            {  
                                for (GlycanSpecies theGlycanSpecies:glycanRxn_compt.getProducts( )){
                                    isProductsinOneCompt= isProductsinOneCompt && theGlycanSpecies.getCompartname().equalsIgnoreCase(
                                            theGlycanCompt.getCompartName());
                                }
                              
                                for (GlycanSpecies theGlycanSpecies:glycanRxn_compt.getReactants()){
                                    isReactantsinOneCompt= isReactantsinOneCompt && theGlycanSpecies.getCompartname().equalsIgnoreCase(
                                            theGlycanCompt.getCompartName());
                                }
                             
                              
                                if(isProductsinOneCompt&&isReactantsinOneCompt){
		    		    theGlycanCompt.getRxns().add(glycanRxn_compt);
		    		     numAdded++;
                                 }
                             }else if(glycanRxn_compt.getType().equalsIgnoreCase("transport")){
                                //influx 
                                if(glycanRxn_compt.getProducts().size()!=0){
                                    for(GlycanSpecies theGlycanSpecies:glycanRxn_compt.getProducts( )){
                                        isProductsinOneCompt = isProductsinOneCompt && theGlycanSpecies.getCompartname().equalsIgnoreCase(
                                            theGlycanCompt.getCompartName());
                                     }
                                    if(isProductsinOneCompt) numInFlux++;
                                } 
                                //outflux  
                                if(glycanRxn_compt.getReactants().size()!=0) {
                                   for (GlycanSpecies theGlycanSpecies: glycanRxn_compt.getReactants()){
                                 /*    if(theGlycanSpecies==null){
                                          System.out.println("size is "+glycanRxn_compt.getReactants().size());
                                          System.out.println("is null");
                                      }
                                  * 
                                  */
                                     isReactantsinOneCompt = isReactantsinOneCompt && theGlycanSpecies.getCompartname().equalsIgnoreCase(
                                           theGlycanCompt.getCompartName());
                                   }
                                    if(isReactantsinOneCompt) numOutFlux++;
                                }

                             }else{
                                numNotAdded++; 
		    	     }
			}
	/*		System.out.println(theGlycanCompt.getRxns().size());
			System.out.println("number of rxns is "+numAdded);
                        System.out.println("number of influx is "+numInFlux);
                        System.out.println("number of outflux is "+numOutFlux);
			System.out.println("number of rxns not added is "+numNotAdded);
	*/		
			
			//add single compartment to vector of compartments
			this.addPathway(theGlycanCompt);
	      }
        }
        
        public void readPathwayInCompts(Model ubModel,boolean reset) 
                throws MalformedURLException, Exception {
		
                int numCompartments =(int) ubModel.getNumCompartments();
		if(reset){
			this.theGlycanRxnNets.clear();
		}			
		
                //create the compartments and add the compartment names
		Vector<GlycanNet> theCompts = new Vector();
                for(int j=0;j<numCompartments;j++){
		    GlycanNet theGlycanCompt = new GlycanNet();
		    theGlycanCompt.setCompartName(this.theCompartments.get(j));
                    theCompts.add(theGlycanCompt);
               }
                
		//assign the species to each compartment
                for(int i = 0; i < (int) this.theGlycanSpecies.size() ; i++)
		{
		    for(GlycanNet theGlycanCompt:theCompts){
                      if(theGlycanSpecies.get(i).getCompartname().equalsIgnoreCase(theGlycanCompt.getCompartName()))
                      {
			theGlycanCompt.addStructure(theGlycanSpecies.get(i));
                        break;
                      }
                    }                    
		}
                
                //add the reactions to each compartment based on compartment information
               for(GlycanNet theGlycanCompt:theCompts){
           //          System.out.println("name of the compartment is "+theGlycanCompt.getCompartName());		
           //          System.out.println("number of the species is "+theGlycanCompt.getStructures().size());
			
		     //add rxns to each single compartment
                        int numInFlux=0,numOutFlux=0,numAdded=0,numNotAdded = 0;
                        
			for(int i=0;i<this.theGlycanRxns.size();i++){
                            boolean isProductsinOneCompt  = true;
                            boolean isReactantsinOneCompt = true; 
			    
                            GlycanRxn glycanRxn_compt = this.theGlycanRxns.get(i);
                                                            
                            if( (glycanRxn_compt.getType().equalsIgnoreCase("reaction")))
                            {  
                                for (GlycanSpecies theGlycanSpecies:glycanRxn_compt.getProducts( )){
                                    isProductsinOneCompt= isProductsinOneCompt && theGlycanSpecies.getCompartname().equalsIgnoreCase(
                                            theGlycanCompt.getCompartName());
                                }
                              
                                for (GlycanSpecies theGlycanSpecies:glycanRxn_compt.getReactants()){
                                    isReactantsinOneCompt= isReactantsinOneCompt && theGlycanSpecies.getCompartname().equalsIgnoreCase(
                                            theGlycanCompt.getCompartName());
                                }
                             
                              
                                if(isProductsinOneCompt&&isReactantsinOneCompt){
		    		    theGlycanCompt.getRxns().add(glycanRxn_compt);
		    		     numAdded++;
                                 }
                             }else if(glycanRxn_compt.getType().equalsIgnoreCase("transport")){
                                //influx 
                                if(glycanRxn_compt.getProducts().size()!=0){
                                    for(GlycanSpecies theGlycanSpecies:glycanRxn_compt.getProducts( )){
                                        isProductsinOneCompt = isProductsinOneCompt && theGlycanSpecies.getCompartname().equalsIgnoreCase(
                                            theGlycanCompt.getCompartName());
                                     }
                                    if(isProductsinOneCompt) numInFlux++;
                                } 
                                //outflux  
                                if(glycanRxn_compt.getReactants().size()!=0) {
                                   for (GlycanSpecies theGlycanSpecies: glycanRxn_compt.getReactants()){
                                 /*    if(theGlycanSpecies==null){
                                          System.out.println("size is "+glycanRxn_compt.getReactants().size());
                                          System.out.println("is null");
                                      }
                                  * 
                                  */
                                     isReactantsinOneCompt = isReactantsinOneCompt && theGlycanSpecies.getCompartname().equalsIgnoreCase(
                                           theGlycanCompt.getCompartName());
                                   }
                                    if(isReactantsinOneCompt) numOutFlux++;
                                }

                             }else{
                                numNotAdded++; 
		    	     }
			}
	/*		System.out.println(theGlycanCompt.getRxns().size());
			System.out.println("number of rxns is "+numAdded);
                        System.out.println("number of influx is "+numInFlux);
                        System.out.println("number of outflux is "+numOutFlux);
			System.out.println("number of rxns not added is "+numNotAdded);
	*/		
			
			//add single compartment to vector of compartments
			this.addPathway(theGlycanCompt);
	      }		
	}

	public void readSBMLfromFile(String sbmlfilename) throws Exception
        {
              try {  
                SBMLReader theSBMLReader = new SBMLReader();
	    	SBMLDocument thesbmlDocument = theSBMLReader.readSBMLFromFile(sbmlfilename);
	        Model theSBMLModel = thesbmlDocument.getModel();
                
	/*      System.out.println("Number of glycan structures in SBML is "+theSBMLModel.getNumSpecies());
	        System.out.println("Number of glycan rxns in SBML is "+theSBMLModel.getNumReactions());
                * 
                */
	        
	        readStructures(theSBMLModel);
                readRxns(theSBMLModel);
                readCompartments(theSBMLModel);                
                readPathwayInCompts();
	    	
	    	int numCompartments =(int) theSBMLModel.getNumCompartments();
	    	for(int i=0;i<numCompartments;i++){
	    		GlycanNet glycoCompt = this.theGlycanRxnNets.get(i);
	   /* 		System.out.println("Compartment "+glycoCompt.getCompartName()+" has "+
	    		                glycoCompt.getStructures().size()+
	    				" species and "+glycoCompt.getRxns().size()+" reactions") ;
                                        * 
                                        */
	    	}	    	
	      }  catch (Exception e) {
	        throw e;
	       }	     
	 }
            
        private GlycanRxn sbml2GlycanRxn(Reaction sbmlRxn,GlycanNet net) throws Exception {
		//	System.out.println("number of species is "+net.getStructures().size());
			
			GlycanRxn theGlycanRxn;
                        GlycanSpecies reac,prod;
			
                        ListOf<SpeciesReference> reactRefList = sbmlRxn.getListOfReactants();
                        ListOf<SpeciesReference> prodRefList = sbmlRxn.getListOfProducts();
                        
                        if(reactRefList.size()> 1 && prodRefList.size()>1){
                           throw new Exception("Current reaction format only supprt one reaction and one product"); 
                        }
                        
                        if(reactRefList.size()==1 && reactRefList.get(0).getSpecies()!= null){ 
                           reac = net.findGlycanSpeciesByName(reactRefList.get(0).getSpecies());
                        } else {
                           reac = null;
                        }
                        
			if(prodRefList.size()==1 && prodRefList.get(0).getSpecies()!= null){
                            prod = net.findGlycanSpeciesByName( prodRefList.get(0).getSpecies());
                        } else {
                            prod = null;
                        }         
                        			
                        if ((reac==null)&&(prod==null)){
				theGlycanRxn = null;
			}else{
				theGlycanRxn = new GlycanRxn(reac,prod);
			}
			return theGlycanRxn;
	}	
	
	private static String note2GlycanStr(String a) {
		a=a.replaceAll("<notes>", "");
		a=a.replaceAll("</notes>", "");
		a=a.replaceAll("<body xmlns=\"http://www.w3.org/1999/xhtml\">", 
                        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
		a=a.replaceAll("</body>", "");
		a=a.trim();
		return a;
	}
        
        private static String annotation2GlycanStr(String a) {
		a=a.replaceAll("<annotation>", "");
		a=a.replaceAll("</annotation>", "");
		a=a.replaceAll("<glycoct xmlns=\"http://www.eurocarbdb.org/recommendations/encoding\">", 
                        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
		a=a.replaceAll("</glycoct>", "");
		a=a.trim();
		return a;
	}
	
	private void addPathways(Collection<GlycanNet> pathways) {
		if ((pathways != null) && (pathways.size() > 0)) {
                    for (GlycanNet pathway : pathways) addPathway(pathway);		 
		}
	 }
         
	private void addPathway(GlycanNet pathway){
           if ((pathway != null) && pathway.getStructures().size()!=0) {
               this.theGlycanRxnNets.add(pathway);
           }
        }
}
	
