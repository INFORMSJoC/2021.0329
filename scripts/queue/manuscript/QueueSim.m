function [MeanTime] = QueueSim(lambda, x, runlength)

    warmup=500;   
    real=36;     
    total=warmup+real;         
    nRuns=runlength;    
    
    interarrival=exprnd(1/x,total,nRuns);          
    service=exprnd(1/lambda,total,nRuns); 
    
    MeanTime=zeros(nRuns,1);
    for nr=1:nRuns                           
        arrival_time = cumsum(interarrival(:,nr));
        leave = zeros(total,1);
        leave(1) = arrival_time(1)+service(1,nr);
        for i=2:total
            leave(i) = max([arrival_time(i),leave(i-1)])+service(i,nr);
        end
        whole = leave - arrival_time;
        MeanTime(nr)=mean(whole((warmup+1):total));     
    end

    
end

