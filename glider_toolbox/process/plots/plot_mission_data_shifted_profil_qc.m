function plot_mission_data_shifted_profil_qc(data,paramplot,figparam)
%  function plot_mission_data_shifted_profil_qc(data,paramplot,figparam)  
%     
% created by L. Houpert (houpertloic@gmail.com), 05/02/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox

     if isfield(data,'qc')
        plotwithqcflag =1;
     else
        plotwithqcflag =0;
     end
     
    disp(['Shifted ' paramplot.plotext ' profiles ' data.gname ] )
     
     xdata   = data.x ;
     ydata   = data.y ;
     xdata00 = data.x00 ;
     ydata00 = data.y00 ;
     time    = data.time;
     
     if plotwithqcflag
        qc_flag = data.qc;
     end
     
     coldef    = figparam.color.coldef;
     nbshiftpf = figparam.shiftplot.shiftpf_nb;

     yaxlim       = paramplot.yaxislim ;
     xref         = paramplot.xref;
     x_shift_val  = paramplot.x_shift_val;
     
     ipf1 = data.ipf1;
     ipf2 = data.ipf2;
     
     
     divephase = data.divephase;

     varshortname = paramplot.varshortname;
     varlongname  = paramplot.varlongname;         
     figrep       = figparam.shiftplot.dir ;
     figsuf       = figparam.shiftplot.name ;
     figext       = paramplot.plotext;
     figtitle     = figparam.shiftplot.titleapp;
     
     fignb = 0; 
     cc=0;
     for ijk = 1:length(ipf1)
        indintref0 = ipf1(ijk):ipf2(ijk);
    	indintref =  indintref0;    
       if length(indintref)>2
          cc = cc +1;
          if mod(cc,nbshiftpf*6)==0 
              savefigure = 1;
          else
              savefigure = 0;
          end        
          if mod(cc,6)==1 % condition to create a new group of profiles
              newgroup = 1;
          else
              newgroup = 0;
          end           
          if mod(cc,6)==0 % condition to to identify the last profile of the group
              lastofgroup = 1;
          else
              lastofgroup = 0;
          end                  
          indint = indintref;%(~isnan(ptmp(indintref)));
          xdatapf = xdata(indint);
          ydatapf = ydata(indint);
          xdatapf00 = xdata00(indint);
          ydatapf00 = ydata00(indint);       
          timepf00  = time(indint);  
          if plotwithqcflag          
            qc_flagpf   = qc_flag(indint);
          end
          divephasepf = divephase(indint);  
          
          idataok   = find(~isnan(xdatapf));        

          if plotwithqcflag      
            inotnan   = find(qc_flagpf~=9); % no nan in the raw data 
          end  
          idatabad  = intersect(find(isnan(xdatapf)),inotnan) ;         
          %---------------------------------------------------------------------------------
          % define the color code in the general case where descent profiles
          % are not always only before or after ascent profiles
          if mod(cc,2) == 1 % first profile should be a descent one (for the color code)
              if median(divephasepf(qc_flagpf~=9))==1 % profile is a dive we don't do anything
                 colundefined = 0;                  
              elseif median(divephasepf(qc_flagpf~=9))== 4 % ascent
                 cc=cc+1; % if after increase index we are in the case where the figure need to be save
                 if savefigure==0 & mod(cc,nbshiftpf*6)==0 
                     savefigure=1;
                 end 
                 if newgroup == 0 & mod(cc,6)==1 
                    newgroup = 1;
                 end      
                 if lastofgroup == 0 &  mod(cc,6)==0
                    lastofgroup = 1;    
                 end
                 colundefined = 0;                    
              else % color code for undefined phase
                 colundefined = 1;
              end
          else % 2nd profile should be a ascent one (for the color code)
              if median(divephasepf(qc_flagpf~=9))==1 % profile is a dive so one profile has to be skipped
                 cc=cc+1;    
                 if savefigure==0 & mod(cc,nbshiftpf*6)==0 
                     savefigure=1;
                 end               
                 if newgroup == 0 & mod(cc,6)==1 
                    newgroup = 1;
                 end
                 if lastofgroup == 0 &  mod(cc,6)==0
                    lastofgroup = 1;    
                 end
                 colundefined = 0;                  
              elseif median(divephasepf(qc_flagpf~=9))== 4 % ascent
                 colundefined = 0;                    
              else % color code for undefined phase
                 colundefined = 1;
              end              
          end   
          % color code for profiles
          if colundefined == 0
              if mod(cc,6)~=0
                col = coldef{mod(cc,6)};
              else
                col = coldef{6};
              end
          else
              col = 'k';
          end
        
        if cc==1
        close all
        fignb = fignb +1 ;
        pl2 = [];
       	figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
       	set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                    'Paperposition',figparam.figpos)
     	ax1=axes('position',figparam.axpos);
    	
        indlegok=[];
        clear pl3ref pl1 pl2 pl3 
	plcolor = [];
        end
           
        set(0,'currentfigure',figg)
        hold on
        
        if lastofgroup % lastgroup actiopn has to appear before the newgroup one
            xtext = shiftx;  
            if strcmp(paramplot.matlabvar_y,'pres') % if profiles are plotted (vertical ax is reversed)
                ytext = maxydata + (yaxlim(2) - yaxlim(1))/10;
                if ytext > yaxlim(2) 
                    ytext  = yaxlim(2) - (yaxlim(2) - yaxlim(1))/10 ;
                end         
                vertalg = 'bottom';           
                ydataref = minydata:ytext;
            else
                ytext = minydata - (yaxlim(2) - yaxlim(1))/10;
                if ytext < yaxlim(1) 
                    ytext  = yaxlim(1) + (yaxlim(2) - yaxlim(1))/10 ;
                end         
                vertalg = 'bottom';           
                ydataref = ytext:(maxydata-ytext)/4:maxydata ;              
            end
            
            pl0=plot(zeros(size(ydataref)) + shiftx ,ydataref,'--','color',figparam.color.greycol{1});   
            
            timeend2nddive = timepf00(end); % time at the end of the group of profiles
            textlabel = [data.gname(1) '-' datestr(timeend2nddive,'dd/mm HH:MM')];
            ht1=text(xtext, ytext , textlabel);
            set(ht1,'horizontalalignment','center','verticalalignment',vertalg)
            set(ht1,'fontsize',figparam.fs);%,'backgroundcolor','w');           
            
        end   
        
        if newgroup 
              ytext = 100; % initialise at the beginning of each new profil    
              xtext = 100;
              % keep track of group property (min, max) for text and or
              % other background lines
              minxdata = min(xdatapf);
              maxxdata = max(xdatapf);
              minydata = min(ydatapf00);
              maxydata = max(ydatapf00);              
       end
        
        
        if min(xdatapf) < minxdata
              minxdata = min(xdatapf);
        end
        if max(xdatapf) > maxxdata;
            maxxdata = max(xdatapf);
        end
         if min(ydatapf00) < minydata
              minydata = min(ydatapf00);
        end
        if max(ydatapf00) >  maxydata;
            maxydata = max(ydatapf00);
        end
                   
          if cc == 1
              shiftx = 0;
              minxdatapf1 = - x_shift_val;     	
              maxxdatapf1 = max((xdatapf00 -xref) + shiftx) + x_shift_val;    
          elseif newgroup & cc>1 
              shiftx = shiftx + 2*x_shift_val;      
          end            

         
        if (max(xdatapf00 -xref) + shiftx) > maxxdatapf1 & length(idataok)>4
            maxxdatapf1 =  (max(xdatapf00 -xref) + shiftx); 
        end
        %if (min(xdatapf00 -xref) + shiftx) < minxdatapf1,minxdatapf1 =  (min(xdatapf00 -xref) + shiftx); end        
                              
        %---- 
                   
        
        ibadqc{1} = idatabad(qc_flagpf(idatabad) == 3); symbqc{1} = 'o'; legstr{1}='flag 3';
        ibadqc{2} = idatabad(qc_flagpf(idatabad) == 4); symbqc{2} = '+'; legstr{2}='flag 4';
        ibadqc{3} = idatabad(qc_flagpf(idatabad) == 8); symbqc{3} = 'd'; legstr{3}='flag 8'; % Seabird basestation flag for interpolated value
        ibadqc{4} = idatabad(qc_flagpf(idatabad) == 41); symbqc{4} = '^'; legstr{4}='flag 41';% spike     
        ibadqc{5} = idatabad(qc_flagpf(idatabad) == 42); symbqc{5} = 'h'; legstr{5}='flag 42'; % gradient test failed
        ibadqc{6} = idatabad(qc_flagpf(idatabad) == 43); symbqc{6} = 'v'; legstr{6}='flag 43'; % density inversion
        ibadqc{7} = idatabad(qc_flagpf(idatabad) == 44); symbqc{7} = '>'; legstr{7}='flag 44'; %  inferior to fix threshold
        ibadqc{8} = idatabad(qc_flagpf(idatabad) == 40); symbqc{8} = 'x'; legstr{8}='flag 40'; %  manually flagged bad     
        ibadqc{9} = idatabad(qc_flagpf(idatabad) == 45); symbqc{9} = '<'; legstr{9}='flag 45'; %  outliers from pressure bin mean statistics           
  
        pl1{cc}=plot((xdatapf00(inotnan) - xref) + shiftx ,ydatapf00(inotnan),'.-','color',figparam.color.greycol{2},'markersize',figparam.mksize2);        
        pl2=plot((xdatapf00(idataok) - xref) + shiftx ,ydatapf00(idataok),'-+','color',col,'markersize',figparam.mksize1);               
        if ~isempty(pl2)
            plcolor =pl2;
        end

        for ipp=1:length(ibadqc)
            
            pl3{ipp,cc} =plot((xdatapf00(ibadqc{ipp}) - xref) + shiftx ,ydatapf00(ibadqc{ipp}),symbqc{ipp},'color',figparam.color.greycol{3},'markersize',figparam.mksize2);      
            if ~isempty(pl3{ipp,cc}) ; 
                pl3ref{ipp}=pl3{ipp,cc};
                indlegok = [indlegok ipp];
            end
        end


        
        if savefigure

              indqfok = unique(indlegok);
              if ~isempty(plcolor) & length(indqfok)>=1
                legcol{1}='flag 0, 1 or 2';
                legend([plcolor pl3ref{indqfok}],[legcol legstr(indqfok)],'location','west')   
              elseif ~isempty(plcolor) & isempty(indqfok)
                legcol{1}='flag 0, 1 or 2';
                legend(plcolor ,legcol ,'location','west')  
              elseif isempty(plcolor) & ~isempty(indqfok)
                legend(pl3ref{indqfok},legstr(indqfok),'location','west')            
              end
              
              
            for izz=1:length(pl1)
                set(pl1{izz},'zdata',get(pl1{izz},'xdata')*0-1)
                for ipp=1:length(legstr)
                    set(pl3{ipp,izz},'zdata',get(pl3{ipp,izz},'xdata')*0+1)            
                end
            end
            
            if strcmp(paramplot.matlabvar_y,'pres') % when y is a pressure variable
                set(gca,'ydir','reverse'); 
            end
            
            
            title([figtitle varlongname ' anomaly (' varshortname ' - ' num2str(xref,'%2.0f') '), shifted by ' num2str(x_shift_val) ''],'fontsize',figparam.fstitle)         
            set(gca,'xlim',[minxdatapf1 maxxdatapf1],...
                'ylim',yaxlim)
            %text(shiftss2(1,3:round(nbshiftpf/6):end),shiftss2(1,3:round(nbshiftpf/6):end)*0-1050,datestr(timedive(iiic2(isel2(3:round(nbshiftpf/6):end))/2),'dd/mm/yy'),'horizontalalignment', 'center')

            set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,             
            print(gcf,'-dpng',[figrep filesep figsuf '_' figext '_' num2str(fignb)]);
            close
   
            % new figure 
            fignb = fignb +1 ;
            pl2 = [];
            figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
            set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                        'Paperposition',figparam.figpos)
            ax1=axes('position',figparam.axpos);

            indlegok=[];
            clear pl3ref pl1 pl2 pl3 
	    plcolor=[];  
            	
            
            shiftx = 0;

            %-----------------------------------------------------------
  
        end
        
      end
     end
end
  
