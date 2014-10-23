package org.glyco;

import java.util.Vector;

public class GlycanRxn {
	private Vector<GlycanSpecies> reactants = new Vector<GlycanSpecies>();
	private Vector<GlycanSpecies> products  = new Vector<GlycanSpecies>();;
	private Enzyme enz;
        private String type;

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }

    /*
     * 
     */    
    public GlycanRxn(GlycanSpecies reac, GlycanSpecies prod) throws Exception
    {
		if(reac!=null)  this.reactants.add(reac);
		if(prod!=null) this.products.add(prod);
                this.type = "reaction";
                if(reac==null||prod==null) {
                    this.type="transport";
                }else if(!reac.getCompartname().equalsIgnoreCase(prod.getCompartname())){
                    this.type="transport";
                }else if(reac.getCompartname().equalsIgnoreCase(prod.getCompartname())){
                    this.type="reaction";
                }else{
                    throw new Exception("Reaction not set up properly");
                }                    
        }
	
    public Vector<GlycanSpecies> getReactants() {
		return reactants;
    }
    
    public void setReactants(Vector<GlycanSpecies> reactants) {
		this.reactants = reactants;
    }
    
    public Vector<GlycanSpecies> getProducts() {
		return products;
    }
	
    public void setProducts(Vector<GlycanSpecies> products) {
		this.products = products;
    }
    
    public Enzyme getEnz() {
		return enz;
    }

    public void setEnz(Enzyme enz) {
		this.enz = enz;
     }
	
     public String toString(){
            
     String reactantsnames="";
     for(int i=0;i<this.reactants.size();i++){
        reactantsnames.concat(this.reactants.get(i).getName());
        if(i!=reactants.size()-1) reactantsnames.concat("_+_");        
     }
     
     String productsnames="";
     for(int i=0;i<this.products.size();i++){
        productsnames.concat(this.products.get(i).getName());
        if(i!=reactants.size()-1) reactantsnames.concat("_+_");        
     } 
     
     if(reactantsnames.isEmpty()&& productsnames.isEmpty()){
	return "unknown name";
     }else{
	return reactantsnames+"_to_"+productsnames;
       }		
     }
}
