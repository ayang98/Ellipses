function [model_ellipse,model_MSE] = RANSAC(x_ran,y_ran)
 
xrange = max(x_ran) - min(x_ran); %x range
yrange = max(y_ran) - min(y_ran); %y range 

%test x data points to see if there are points outside the normal range of an ellipse 
if xrange > 75 && min(x_ran) < 50 %if we have a number of points far from ellipse 
   include = setdiff(1:length(x_ran),find(x_ran < 50)); %only include points closer to the ellipse
   x_ran = x_ran(include);
   y_ran = y_ran(include);
   
   xrange = max(x_ran) - min(x_ran); %redefine x range for next if statement 
   yrange = max(y_ran) - min(y_ran); %redefine y range for next set of if statements
end

if xrange > 75 && max(x_ran) > 140 
   include = setdiff(1:length(x_ran),find(x_ran > 140)); 
   x_ran = x_ran(include);
   y_ran = y_ran(include);
   
   xrange = max(x_ran) - min(x_ran);
   yrange = max(y_ran) - min(y_ran); %redefine y range for next set of if statements
end

%adjust points if set contains large range of y points, signaling points
%outside of an ellipse
if yrange > 75 && min(y_ran) < 120
    include = setdiff(1:length(y_ran),find(y_ran < 120)); 
    x_ran = x_ran(include);
    y_ran = y_ran(include);
    yrange = max(y_ran) - min(y_ran); %redefine y range for next if statement
end

if yrange > 75 && max(y_ran) > 170
   include = setdiff(1:length(y_ran),find(y_ran > 170));
   x_ran = x_ran(include);
   y_ran = y_ran(include);
end

%fit initial ellipse
[init_model] = fit_ellipse(x_ran,y_ran);

%axes of original fitted ellipse
xr = init_model.a;
yr = init_model.b;
%ecc = axes2ecc(max(xr,yr),min(xr,yr)); %eccentricity
ecc = 0.4;
%standard deviations
x_stand_dev = std(x_ran); 
y_stand_dev = std(y_ran);

%averages of the points
x_av = mean(x_ran(:));
y_av = mean(y_ran(:));

if ecc > 0.75 && yr > xr %if the ellipse is stretched horizontally by a large extent
   x_low = find(x_ran < x_av - 1.25*x_stand_dev); %find low points that should be excluded
   x_high = find(x_ran > x_av + 1.25*x_stand_dev); %find high x points that should be excluded
   exclude = union(x_low,x_high); %combine points that should be excluded
   include = setdiff(1:length(x_ran),exclude); %find points that should be included
   y_ran = y_ran(include); %included x
   x_ran = x_ran(include); %included y
elseif ecc > 0.75 && xr > yr  %if the ellipse is stretched vertically by a large extent
   y_low = find(y_ran < y_av - 1.25*y_stand_dev); %find low y points that should be excluded
   y_high = find(y_ran > y_av + 1.25*y_stand_dev); %find high y points that should be excluded
   exclude = union(y_low,y_high); %combine excluded points
   include = setdiff(1:length(y_ran),exclude);
   x_ran = x_ran(include); %included x
   y_ran = y_ran(include); %included y
end


%%%RANSAC 

%Variables user must currently define threshold 
nis = 7; %number of values in each sample
maxiter = 600;

best_error = +Inf; %initialize best error for our created models
tolerance = 5; %error tolerance for point from model 
max_con_set = 0.7*length(x_ran); %minimum requirement for parameter fitting

[model_ellipse] = fit_ellipse(x_ran,y_ran); %initial model in case others don't w

for i = 1:maxiter %for each iteration 
    
    vals = randi(length(x_ran),1,nis); %randomly determine index values for initial test
    samplex = x_ran(vals); %x entries from randomly selected indices
    sampley = y_ran(vals); %y 

    test_ind = setdiff(1:length(x_ran),vals);   
    Xtest = x_ran(test_ind); %entries that are in x_ran but not x_set for testing against model
    Ytest = y_ran(test_ind);

    consset = []; %initialize consensus set for the model 
    
    [sample_model] = fit_ellipse(samplex,sampley); %create model structure from sample
    if ~isempty(sample_model)
    [~,x_model,y_model] = ellipse(sample_model.a,sample_model.b,-sample_model.phi,sample_model.X0_in,sample_model.Y0_in,0,'k',[]); %return test coordinates from structure  
    end

    %the model for this sample has now been defined
    
    for j = 1:length(Xtest) %for each point of test points
        
        mindist = +Inf; %initialize minimum distance for test point from model as infinite 
       
        for k = 1:length(x_model) %for each point of model points    
            dist = sqrt((Xtest(j)-x_model(k))^2 + (Ytest(j)-y_model(k))^2); %calculate distance between test and model point                       
            if dist < mindist %if new calculated distance is smaller than current minimum distance from ellipse
               mindist = dist; %this is the new minimum distance
            end            
        end %end after testing all points against points of model ellipse
        
        if mindist < tolerance %if within threshold
            consset_new = [Xtest(j); Ytest(j)]; 
            consset = [consset_new consset]; %add this point to consensus set 
        end
    end
    
    if size(consset,2) >= max_con_set %if we have enough values in consensus set, we have a reasonably good model
       [new_model] = fit_ellipse(consset(1,:),consset(2,:)); %create model structure from consensus set
       [~,x_model,y_model] = ellipse(new_model.a,new_model.b,-new_model.phi,new_model.X0_in,new_model.Y0_in,0,'k',[]); %reestimate model with parameters
       
        model_MSE = 0; %initialize model error for this model
        for j = 1:size(consset,2) %for each point in consensus set 
            mindist = +Inf; %initialize minimum distance as large number
            for k = 1:length(x_model) %for each point in model
                dist = sqrt((consset(1,j) - x_model(k))^2 + (consset(2,j) - y_model(k))^2); %calculate distance from 
                if dist < mindist %if newly calculated distance is smaller than current minimum
                   mindist = dist; %this is the new minimum for that point
                end
            end
            model_MSE = model_MSE + mindist^2; %add distance to mean-squared error
        end 
        model_MSE = model_MSE/size(consset,2); %calculate the mean part of mean-squared error
        
        if model_MSE < best_error %if the new model fits the consensus set points better
           model_ellipse = new_model; %accept as a better model
           best_error = model_MSE; %we have a new best error 
        end %ending model acceptance/rejection
            
    end %end if statement as to whether or not the model of this iteration is reasonably good
    
end %all models necessary have been created

model_MSE = roundn(best_error,0); %round error to nearest whole number
end