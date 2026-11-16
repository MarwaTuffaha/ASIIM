classdef plotProp
    %plotProp Plot options for figures genorated from a assim object
    %   Detailed explanation goes here

    properties
        fileName
        region = '';
        ageRange = 1:16;
        visible = 1;
        showTitle = 0;
        save = 0;
        subPlot = 0;
        percent = 1;
        fontSize = 18;
        lineWidth = 1.3;
        verticalLinesFlag = true;
        plotSize = 1200;
        I4Only = false;
        showLegend = false;
        yLabel = 'y';
        yLim = [-inf,inf];
    end
    properties (Hidden)
        goldenRatio = ( 1 + sqrt(5) ) / 2;
    end
    properties (Dependent)
        position
    end

    methods
        function position = get.position(obj)
            position = [10 40 obj.plotSize*obj.goldenRatio obj.plotSize];
        end
    end
end

