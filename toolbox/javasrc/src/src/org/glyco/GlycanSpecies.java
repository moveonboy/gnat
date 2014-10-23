package org.glyco;

import java.io.Serializable;
import java.util.HashMap;
import org.eurocarbdb.application.glycanbuilder.Glycan;

public class GlycanSpecies implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
        protected boolean isGlycan;
	protected Glycan  structure;
	protected String  name;
        protected String id;
	private String compartname;
	private HashMap<String,String> dbIDs;
	
	public Glycan getStructure() {
		return structure;
	}
	public void setStructure(Glycan structure) {
		this.structure = structure;
	}
	public String getName(){
		return name;
	}
	public void setName(String name){
		this.name = name;
	}
	
	public GlycanSpecies(Glycan struct){
              this(struct,struct.toString());
        }
                
        public GlycanSpecies(Glycan struct,String name){
		this.structure = struct;
		this.name = name;
		this.setDbIDs(new HashMap<String, String>());
	}
        
        public GlycanSpecies(Glycan struct,String name,String comptName){
		this.structure = struct;
		this.name = name;
                this.compartname = comptName;
		this.setDbIDs(new HashMap<String, String>());
	}        
        
	public String getCompartname() {
		return compartname;
	}
	public void setCompartname(String compartname) {
		this.compartname = compartname;
	}
	public HashMap<String,String> getDbIDs() {
		return dbIDs;
	}
	public void setDbIDs(HashMap<String,String> dbIDs) {
		this.dbIDs = dbIDs;
	}
	
	public void addDbID(String resourceName,String resourceID){
		this.dbIDs.put(resourceName,resourceID);
	}

    public String getID() {
        return id;
    }
    
    public void setID(String id){
        this.id = id;
    }

}
