function v1 = ab_method(v0,v1,h,func)


    
v1 = ((h/2)*((3*func(v1))-func(v0)))+v1;

end
    

