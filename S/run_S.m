for runCount=1:3
    mkdir('H:\Documents\MATLAB\Results\S\Lin',strcat('Data',num2str(runCount)));
    mkdir('H:\Documents\MATLAB\Results\S\Lin',strcat('Img',num2str(runCount)));
    cd(strcat('H:\Documents\MATLAB\Results\S\Lin\Data',num2str(runCount)));
    run learnImp_S_Lin.m
    
    cd(strcat('H:\Documents\MATLAB\Results\S\Lin\Img',num2str(runCount)));    
    run saveImages.m;
    
    %%
    mkdir('H:\Documents\MATLAB\Results\S\GP10',strcat('Data',num2str(runCount)));
    mkdir('H:\Documents\MATLAB\Results\S\GP10',strcat('Img',num2str(runCount)));
    cd(strcat('H:\Documents\MATLAB\Results\S\GP10\Data',num2str(runCount)));
    run learnImp_S_GP10.m
    
    cd(strcat('H:\Documents\MATLAB\Results\S\GP10\Img',num2str(runCount)));    
    run saveImages.m;
    
    %%
    mkdir('H:\Documents\MATLAB\Results\S\GP10',strcat('Data',num2str(runCount)));
    mkdir('H:\Documents\MATLAB\Results\S\GP10',strcat('Img',num2str(runCount)));
    cd(strcat('H:\Documents\MATLAB\Results\S\GP25\Data',num2str(runCount)));
    run learnImp_S_GP25.m
    
    cd(strcat('H:\Documents\MATLAB\Results\S\GP25\Img',num2str(runCount)));    
    run saveImages.m;
    
    %%
    mkdir('H:\Documents\MATLAB\Results\S\Lin',strcat('Data',num2str(runCount)));
    mkdir('H:\Documents\MATLAB\Results\S\Lin',strcat('Img',num2str(runCount)));
    cd(strcat('H:\Documents\MATLAB\Results\S\Lin\Data',num2str(runCount)));    
    run learnImp_S_simple.m
    
    cd(strcat('H:\Documents\MATLAB\Results\S\Lin\Img',num2str(runCount)));   
    run saveImages.m;
    
    %%
    clc;
    fprintf('made it');
end
