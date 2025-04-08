function [SigStaticClean,TtestvGvalClean] = CleanSigMatrix(SigStatic,testvGval,NumConsecElim,clusters)

Tbin = numel(SigStatic(:,1));
bconseq = 0;
SigStaticClean =  SigStatic;
testvGvalClean = testvGval;


for bin = 1:Tbin    
    for bin2 = 1:Tbin-1       
        if bin == bin2 %removinf diagonal
            %             if bin2 > 3
            %                 lo = 3;
            %             else
            %                 lo = bin2-1;
            %             end
            %             if bin2 > Tbin-3
            %                 hi = 1;
            %             else
            %                 hi = 2;
            %             end
            %             SigStaticClean(bin,bin2-lo:bin2+hi) = 0;
            %             SigStatic(bin,bin2-lo:bin2+hi) = 0;
            %             for ii = bin2-lo:bin2+hi
            %                testvGvalClean{bin,ii} = NaN;
            %             end
             SigStaticClean(bin,bin2) = 0;
       %%     SigStatic(bin,bin2) = 0;
       %%     testvGvalClean{bin,bin2} = NaN;
        end
    end
end
ii = 61;
for bin = 1:Tbin    
    for bin2 = 1:Tbin
       TtestvGvalClean(bin,bin2) = nanmean(testvGvalClean{bin,bin2});
       TtestvGvalClean2(ii,bin2) = nanmean(testvGvalClean{bin,bin2});       
    end
    ii = ii-1;
end


for bin = 1:Tbin-1
    i=1;
    bconseq = 0;
    for bin2 = 2:Tbin-1
        
             
        bsig = SigStatic(bin,bin2);
        bsig2 = SigStatic(bin,bin2+1);
                
        if bsig & bsig2 == 0 %eliminating 1 sig bin
            bincon(i) = bin2;
            maxcon(i) = bconseq+1; 
            bconseq = 0;
            i = i +1;
        end
        if bsig & bsig2
            bconseq = bconseq+1;           
        else
            if bconseq
                bincon(i) = bin2;
                maxcon(i) = bconseq+1;
                bconseq = 0;
                i = i +1;
            end
        end
    end
    if i > 1
        for ii = 1:i-1
            if maxcon(ii) < NumConsecElim
                bin
                SigStaticClean(bin,bincon(ii)-maxcon(ii):bincon(ii)) = 0;
               
            end
        end
    end
 
    clear bincon maxcon ii i
    
end


figure
hold on 
imagesc(SigStaticClean);
xlabel('Train bin #','color', 'k')
ylabel('Test bin #','color', 'k')
Z = zeros(Tbin+1);
hm = mesh(0.5:1:61.5,0.5:1:61.5,Z);
hm.FaceColor = 'none';
hm.EdgeColor = 'k';
xp = 1:61;
for ii = 1:61
    yp(ii) = 31;
end
plot(xp,yp, '-w', 'linewidth',1.3);

yp = 1:61;
for ii = 1:61
    xp(ii) = 31;
end
plot(xp,yp, '-w', 'linewidth',1.3);
axis([0.5,61.5,0.5,61.5])
%axis([1,61,1,61])
set(gca, 'TickDir', 'out');                     %  switch side of axis for tick marks
set(gca, 'TickLength', [.01 .01],'linewidth',1.3);            %  this in proportion of x axis%
set(gca, 'XTick', [1:10:61]);        %  reset tick spacing
set(gca, 'XTickLabel', [1:10:61]);        %  reset tick spacing
set(gca, 'YTick', [1:10:61]);        %  reset tick spacing
caxis([0 1])
axis square; %square
hold off



 
 % cmap = colormap;
figure
colormap('jet');
imagesc(TtestvGvalClean2);
xlabel('Train bin #','color', 'k')
ylabel('Test bin #','color', 'k')
caxis([0 1])
axis square; %square



