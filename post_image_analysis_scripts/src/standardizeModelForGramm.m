function model = standardizeModelForGramm(model)

%     modelStr = func2str(model);
%     modelStr = strrep(modelStr, 'd', 'x');
%     modelStr = strrep(modelStr, 'p(1)', 'a');
%     modelStr = strrep(modelStr, 'p(2)', 'b');
%     modelStr = strrep(modelStr, 'p(3)', 'c');
%     modelStr = strrep(modelStr, 'p(4)', 'd');
%     modelStr = strrep(modelStr, '@(p,', '@(a,b,c,d,');
%     model = str2func(modelStr);
    
    modelStr = func2str(model);
    modelStr = strrep(modelStr, 'd', 'x');
    modelStr = strrep(modelStr, 'p(1)', 'amplitude');
    modelStr = strrep(modelStr, 'p(2)', 'KD');
    modelStr = strrep(modelStr, 'p(4)', 'offset');

        if ~contains(modelStr, 'p(5)')
            modelStr = strrep(modelStr, 'p(3)', 'n');
            modelStr = strrep(modelStr, '@(p,', '@(amplitude,KD,n,offset,');
        elseif contains(modelStr, 'p(5)')
            modelStr = strrep(modelStr, 'p(3)', 'n'); %[pol]/KD_pol
            modelStr = strrep(modelStr, 'p(5)', 'omegaDP'); %pol-dorsal interaction term
            modelStr = strrep(modelStr, '@(p,', '@(amplitude,KD,n,offset, omegaDP,');
        end
        model = str2func(modelStr);

end