package org.glyco;

import java.util.Vector;

public class GlycanNet {
	protected Vector<GlycanSpecies> structures;
	protected Vector<GlycanRxn> rxns;
	private String compartName;
		
	public GlycanNet(){
		this.structures= new Vector<GlycanSpecies>();
		this.rxns = new Vector<GlycanRxn>();
	}
	
	public GlycanNet(String compartName) {
		super();
		this.compartName = compartName;
	}

	public Vector<GlycanSpecies> getStructures() {
		return structures;
	}
	public void setStructures(Vector<GlycanSpecies> structures) {
		this.structures = structures;
	}
	public Vector<GlycanRxn> getRxns() {
		return rxns;
	}
	public void setRxns(Vector<GlycanRxn> rxns) {
		this.rxns = rxns;
	}
	
	public boolean addStructure(GlycanSpecies structure){
		
		if(!this.structures.contains(structure)){
			 this.structures.add(structure);
			 return true;
	     }else{
	    	 return false;
	     }	
	}
	
	
	
	public void addStructures(Vector<GlycanSpecies> structures){
		for(int i=0;i<structures.size();i++){
		  addStructure(structures.get(i));  	
		}		
	}	
	
	public boolean addRxn(GlycanRxn rxnToAdd){
	    if(rxnToAdd==null) return false;
		
		if(!this.rxns.contains(rxnToAdd)){
		    
                    this.rxns.add(rxnToAdd);
		    
                    if(!this.structures.containsAll(rxnToAdd.getReactants())){
		    	for(GlycanSpecies theGlycan: rxnToAdd.getReactants()){
                           if(!this.structures.contains(theGlycan)) this.structures.add(theGlycan);
                        }
		    }		
		    
                    if(!this.structures.contains(rxnToAdd.getProducts())){
                        for(GlycanSpecies theGlycan: rxnToAdd.getProducts()){
                           if(!this.structures.contains(theGlycan)) this.structures.add(theGlycan);
                        }
                    } 
                      
		    return true;
		} else {
		    return false;
		}			
	}
	
	public boolean addRxn(GlycanSpecies reac,GlycanSpecies prod) throws Exception{
		//assuming glycan will be added to structures first
		GlycanRxn rxnToAdd = new GlycanRxn(reac,prod); 
		return this.addRxn(rxnToAdd);
	}
	
	public void addRxns(Vector<GlycanRxn> rxnsToAdd){
		for(int i=0;i<rxnsToAdd.size();i++){
		  if(!this.rxns.contains(rxnsToAdd.get(i))) addRxn(rxnsToAdd.get(i)); 
		}
	}	
	
	public GlycanSpecies findGlycanSpeciesByName(String glycanName){
		for(int i=0;i<this.structures.size();i++){
		    if(this.structures.get(i).getName().equals(glycanName))	return this.structures.get(i);		    
		}
		return null;
	}

	public String getCompartName() {
		return compartName;
	}

	public void setCompartName(String compartName) {
		this.compartName = compartName;
	}
	
	public int findMaxNumResidues(){
		int maxNumResidues=0;
		int numResidues;
		for(int i=0;i<this.structures.size();i++){
			numResidues=this.structures.get(i).getStructure().getCount();
			if(numResidues>maxNumResidues) maxNumResidues=numResidues;
		}
		return maxNumResidues;		
	}
	
	public int findMinNumResidues(){
		int minNumResidues=1000;
		int numResidues;
		for(int i=0;i<this.structures.size();i++){
			numResidues=this.structures.get(i).getStructure().getCount();
			if(numResidues<minNumResidues) minNumResidues=numResidues;
		}
		return minNumResidues;		
	}
	
	public Vector<GlycanSpecies> findGlycanSpeciesByResidueNum(int numResToFind){
		Vector<GlycanSpecies> results=new Vector<GlycanSpecies>();
		int numResidues;
		for(int i=0;i<this.structures.size();i++){
			numResidues=this.structures.get(i).getStructure().getCount();
			if(numResidues==numResToFind) results.add(this.structures.get(i));
		} 
		return results;		
	}
	
	public int findmaxNumElements(){
		int minNumResidues = this.findMinNumResidues();
		int maxNumResidues = this.findMaxNumResidues();
		
		int maxNumElements = 0;
		for(int i=minNumResidues;i<maxNumResidues+1;i++){
			Vector<GlycanSpecies> glycanArrayByResidueNum = findGlycanSpeciesByResidueNum(i);
			if(maxNumElements<glycanArrayByResidueNum.size()) maxNumElements=glycanArrayByResidueNum.size();
		}
		return maxNumElements;
	}

	public GlycanNet clone(boolean b) {
		return null;
	}	
}
