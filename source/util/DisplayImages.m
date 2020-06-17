% ---------------------------------------------------------------------------
% **********************CGRI SIMULATION TOOLBOX***************************
% File: DisplayImages.m
% Purpose: A matlab utility to display an array of images
% Inputs:
% Notes: 
%          - gridsize = [gridx,gridy] is the display grid size
%          - u is an array of size [nx,ny,nz] where each slice is an image
%          - options: DisplayImages(...,'stride',s,...) : plot every s-th image
%                     DisplayImages(...,'start',a,...) : start w/ a-th
%                     image
%                     DisplayImages(...,'gap',g,...): put a gap of size g
%                       in between every image (normalized units e.g. 0.01)
%                     DisplayImages(...,'x',xlabel,'y',ylabel): use the
%                       vectors xlabel and ylabel to label the x and y axes.
%                     DisplayImages(...,'dispaxis',true) displays the axes
%                     DisplayImages(...,'axisunits','cm') displays a string
%                       for the axis units.
%                     DisplayImages(...,'colorbar',true) displays a
%                       colorbar for *each* image (need to add ability to
%                       have a *single* colorbar...)
%                     DisplayImages(...,'margin',[ml,mr,mt,mb]) gives the
%                       images a margin to allow for additional text
%                     DisplayImages(...,'normalize',true) gives all the
%                       plots the same color scale [cmin,cmax]
%                     DisplayImages(...,'clim',[cmin,cmax]) gives every plot
%                       the color scale [cmin,cmax]
% Author:  Nick Henscheid
% Date:    9-2016
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without technical
% support, and with no warranty, express or implied, as to its usefulness for
% any purpose.
% ---------------------------------------------------------------------------

function fig = DisplayImages(gridsize,u,varargin)
    imgn   = size(u,1);
    imgm   = size(u,2);
    imgk   = size(u,3);

    gridy  = gridsize(1);
    gridx  = gridsize(2);
    gridtotal = gridx*gridy;
    aspect = gridx/gridy;

    scrsz = get(groot,'ScreenSize');

    default_gap = 5;  % pixels
    default_stride = 1;checkstride = @(x)(x==floor(x));
    default_x = linspace(0,1,10);
    default_y = linspace(0,1,10);
    p = inputParser;
    addParameter(p,'gap',default_gap,@isnumeric);
    addParameter(p,'start',1,@isnumeric);
    addParameter(p,'stride',default_stride,checkstride);
    addParameter(p,'x',default_x,@isnumeric);
    addParameter(p,'y',default_y,@isnumeric);
    addParameter(p,'axis',false,@islogical);
    addParameter(p,'axisunits','',@ischar)
    addParameter(p,'colorbar',false,@islogical);
    addParameter(p,'margin',[0,0,0,0],@isnumeric);
    addParameter(p,'normalize',false,@islogical); % Makes all color scales the same
    addParameter(p,'clim',[0,1],@isnumeric); % Sets custom common color range
    addParameter(p,'crel',false,@islogical); % Sets color relative to each image
    addParameter(p,'title',{},@ischar); % Cell array of titles to use
    parse(p,varargin{:});
    gap = p.Results.gap;
    start = p.Results.start;
    stride = p.Results.stride;
    x = p.Results.x;
    y = p.Results.y;
    dispaxis = p.Results.axis;
    axisunits = p.Results.axisunits;
    cbar = p.Results.colorbar;
    normalize = p.Results.normalize;
    T = p.Results.title;
    clim = p.Results.clim;
    crel = p.Results.crel;
    cmin = clim(1);
    cmax = clim(2);
    if(normalize)
        idx = start:stride:gridtotal;
        cmin = min(min(min(u(:,:,idx))));
        cmax = max(max(max(u(:,:,idx))));
    end
    marg = p.Results.margin;
    wextra = marg(1)+marg(2);
    hextra = marg(3)+marg(4);
    
    subimgx = (1-(gridx+1)*gap)/gridx;
    subimgy = (1-(gridy+1)*gap)/gridy;

    if(gridx>=gridy)
        % Figure is short and fat or square; set width first.
        figw  = 0.8*scrsz(4);
        figh  = figw/aspect;
        figl  = scrsz(3)/2-figw/2;
        figb  = scrsz(4)/2-figh/2;
        subimgx = (figw-(gridx+1)*gap)/gridx;
        subimgy = (figh-(gridy+1)*gap)/gridy;
        fig=figure('Position',[figl,figb,figw+wextra,figh+hextra]);
    else
        % Figure is tall and skinny; set height first.
        figh  = 0.8*scrsz(4);
        figw  = aspect*figh;
        figl  = scrsz(3)/2-figw/2;
        figb  = scrsz(4)/2-figh/2;
        subimgx = (figw-(gridx+1)*gap)/gridx;
        subimgy = (figh-(gridy+1)*gap)/gridy;
        fig=figure('Position',[figl,figb,figw+wextra,figh+hextra]);
    end
    
    idx = start;
    idt = 1;
    subaxes = cell(gridx,gridy);
    for iy = 1:gridy
        bot = (gridy-iy)*subimgy + (gridy-iy+1)*gap + marg(3);
        for ix = 1:gridx
            left = (ix-1)*subimgx + ix*gap + marg(1);
            subaxes{ix,iy} = axes('units','pixels','position',[left,bot,subimgx,subimgy-30]);
            if(idx<=imgk)
                if(crel)
                    cmin = min(min(u(:,:,idx)))
                    cmax = max(max(u(:,:,idx)))
                    imagesc(x,y,u(:,:,idx),[cmin,cmax]);
                    axis image;
                    set(gca,'YDir','normal');
                else
                    imagesc(x,y,u(:,:,idx),[cmin,cmax]);
                    axis image;
                    set(gca,'YDir','normal');
                end
                
                xtick = linspace(min(x),max(x),5)';
                ytick = linspace(min(x),max(x),5)';
                xlabels = cellstr([num2str(xtick,'%1.2f '),repmat(axisunits,[5,1])]);
                ylabels = cellstr([num2str(ytick,'%1.2f '),repmat(axisunits,[5,1])]);
                set(gca,'XTick',xtick,'YTick',ytick,'XTickLabel',xlabels,'YTickLabel',ylabels,'FontSize',16);
                if(~dispaxis)
                    axis off;
                end 
                if(cbar)
                    colorbar;
                end
                % Old title method - didn't leave myself enough space, fix
                % later
%                 if(numel(T)>0)
%                     title(subaxes{ix,iy},T{idt},'interpreter','latex','fontsize',16);
%                 end
            end
            idx = idx + stride;
            idt = idt + 1;
        end
    end
    
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.98,T,'FontSize',20,'HorizontalAlignment','center');


end