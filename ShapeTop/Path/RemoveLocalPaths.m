function [] = RemoveLocalPaths()
rmpath('./Domains');
rmpath('./Parse');
rmpath('./Fem');
rmpath('./Opt');
rmpath('./Plot');
end
