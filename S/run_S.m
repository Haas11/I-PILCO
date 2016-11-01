for runCount=1:3
    %%
    mkdir('H:\Documents\MATLAB\Results\S\GP10',strcat('Data',num2str(runCount)));
    mkdir('H:\Documents\MATLAB\Results\S\GP10',strcat('Img',num2str(runCount)));
    cd(strcat('H:\Documents\MATLAB\Results\S\GP10\Data',num2str(runCount)));
    run learnImp_S_GP10.m
    
    cd(strcat('H:\Documents\MATLAB\Results\S\GP10\Img',num2str(runCount)));
    run saveImages.m;
    
    %%
<<<<<<< HEAD
    runcount=4;
%     mkdir('H:\Documents\MATLAB\Results\S\GP10',strcat('Data',num2str(runCount)));
%     mkdir('H:\Documents\MATLAB\Results\S\GP10',strcat('Img',num2str(runCount)));
%     cd(strcat('H:\Documents\MATLAB\Results\S\GP25\Data',num2str(runCount)));
    run learnImp_S_GP25.m
    
%     cd(strcat('H:\Documents\MATLAB\Results\S\GP25\Img',num2str(runCount)));    
=======
    mkdir('H:\Documents\MATLAB\Results\S\GP25',strcat('Data',num2str(runCount)));
    mkdir('H:\Documents\MATLAB\Results\S\GP25',strcat('Img',num2str(runCount)));
    cd(strcat('H:\Documents\MATLAB\Results\S\GP25\Data',num2str(runCount)));
    run learnImp_S_GP25.m
    
    cd(strcat('H:\Documents\MATLAB\Results\S\GP25\Img',num2str(runCount)));
>>>>>>> d5abe78d437e0f83a378b815e55750f3b2c4c9c6
    run saveImages.m;
    runcount=4;
     
    %%
    mkdir('H:\Documents\MATLAB\Results\S\simple',strcat('Data',num2str(runCount)));
    mkdir('H:\Documents\MATLAB\Results\S\simple',strcat('Img',num2str(runCount)));
    cd(strcat('H:\Documents\MATLAB\Results\S\simple\Data',num2str(runCount)));
    run learnImp_S_simple.m
    
    cd(strcat('H:\Documents\MATLAB\Results\S\Lin\Img',num2str(runCount)));
    run saveImages.m;
    %%
    mkdir('H:\Documents\MATLAB\Results\S\Lin',strcat('Data',num2str(runCount)));
    mkdir('H:\Documents\MATLAB\Results\S\Lin',strcat('Img',num2str(runCount)));
    cd(strcat('H:\Documents\MATLAB\Results\S\Lin\Data',num2str(runCount)));
    run learnImp_S_Lin.m
    
    cd(strcat('H:\Documents\MATLAB\Results\S\Lin\Img',num2str(runCount)));
    run saveImages.m;
    
    %%
    clc;
    fprintf('made it');
end
