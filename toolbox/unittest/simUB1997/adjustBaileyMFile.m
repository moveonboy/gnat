function adjustMFile(mfilename, parametername, parameternevalue) 
    mfilename = [mfilename '.m'];
    mfilestr  = fileread(mfilename);
    
    %change q value
    strrep   = [parametername '\w*\s*=\s\w*(\.\w*)?;' ];
    obj      = regexp(mfilestr,strrep,'match');
    for i=1:length(obj)
    str = obj{i};
    strberep  = regexp(str,'=\s\w*(\.\w*)?;', 'match');
    strtorep  = ['= ' num2str(parameternevalue) ';'];
    objstr    = regexprep(str,strberep,strtorep);
    mfilestr=regexprep(mfilestr,str,objstr);
    end
    fid = fopen(mfilename,'w');
    fprintf(fid,'%s',mfilestr);
    fclose(fid);
   rehash toolboxcache ;
end