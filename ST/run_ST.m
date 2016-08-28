for runCount=1:3   
    %%
    mkdir('D:\victor\MATLAB\Results\ST\GP25',strcat('Data',num2str(runCount)));
    mkdir('D:\victor\MATLAB\Results\ST\GP25',strcat('Img',num2str(runCount)));
    cd(strcat('D:\victor\MATLAB\Results\ST\GP25\Data',num2str(runCount)));
    run learnImp_ST_GP25.m
    
    cd(strcat('D:\victor\MATLAB\Results\ST\GP25\Img',num2str(runCount)));    
    run saveImages.m;
    
        %%
    mkdir('D:\victor\MATLAB\Results\ST\GP50',strcat('Data',num2str(runCount)));
    mkdir('D:\victor\MATLAB\Results\ST\GP50',strcat('Img',num2str(runCount)));
    cd(strcat('D:\victor\MATLAB\Results\ST\GP50\Data',num2str(runCount)));
    run learnImp_ST_GP50.m
    
    cd(strcat('D:\victor\MATLAB\Results\ST\GP50\Img',num2str(runCount)));    
    run saveImages.m;
    %%
    mkdir('D:\victor\MATLAB\Results\ST\GP50',strcat('DataK2',num2str(runCount)));
    mkdir('D:\victor\MATLAB\Results\ST\GP50',strcat('ImgK2',num2str(runCount)));
    cd(strcat('D:\victor\MATLAB\Results\ST\GP50\Data',num2str(runCount)));    
    run learnImp_STk2_GP50.m
    
    cd(strcat('D:\victor\MATLAB\Results\ST\GP50\ImgK2',num2str(runCount)));   
    run saveImages.m;
    
    %%
    clc;
    fprintf('made it');
end
