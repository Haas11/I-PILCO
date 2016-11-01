for runCount=1:3   
    %%
    mkdir('D:\victor\MATLAB\Results\ST\GP75',strcat('DataK3_',num2str(runCount)));
    mkdir('D:\victor\MATLAB\Results\ST\GP75',strcat('ImgK3_',num2str(runCount)));
    cd(strcat('D:\victor\MATLAB\Results\ST\GP75\DataK3_',num2str(runCount)));    
    run learnImp_STk2_GP50.m
    
    cd(strcat('D:\victor\MATLAB\Results\ST\GP75\ImgK3_',num2str(runCount)));   
    run saveImages.m;
    
    %%
    runCount = 4;
    mkdir('D:\victor\MATLAB\Results\ST\GP25',strcat('Data',num2str(runCount)));
    mkdir('D:\victor\MATLAB\Results\ST\GP25',strcat('Img',num2str(runCount)));
    cd(strcat('D:\victor\MATLAB\Results\ST\GP25\Data',num2str(runCount)));
    run learnImp_ST_GP25.m 
    
    cd(strcat('D:\victor\MATLAB\Results\ST\GP25\Img',num2str(runCount)));    
    run saveImages.m;
    runCount = 4;

        %%
    runCount = 0
    mkdir('D:\victor\MATLAB\Results\ST\GP50',strcat('Data',num2str(runCount)));
    mkdir('D:\victor\MATLAB\Results\ST\GP50',strcat('Img',num2str(runCount)));
    cd(strcat('D:\victor\MATLAB\Results\ST\GP50\Data',num2str(runCount)));
    run learnImp_ST_GP50.m
    
    cd(strcat('D:\victor\MATLAB\Results\ST\GP50\Img',num2str(runCount)));    
    run saveImages.m;

    runCount = 0;
    
    %%
    clc;
    fprintf('made it');
end
