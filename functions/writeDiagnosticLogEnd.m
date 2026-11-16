function writeDiagnosticLogEnd(diagnosticFileID, f, initKval, objStart, thisFitTiming, fitTime)
    % Log the objective function value at the final point
    objEnd = f(initKval);
    objImprovement = objEnd / objStart;
    
    fprintf(diagnosticFileID, '----------------------\n');
    fprintf(diagnosticFileID, 'Objective Function Value at Final Point: %f\n', objEnd);
    fprintf(diagnosticFileID, 'Objective Function Improvement Ratio: %f\n', objImprovement);
    fprintf(diagnosticFileID, 'Total Fit Timing for This Step: %f\n', thisFitTiming);
    fprintf(diagnosticFileID, 'Cumulative Fit Timing: %f\n', fitTime);
    fprintf(diagnosticFileID, 'Diagnostic log ended on %s\n', datestr(now));
    
    fclose(diagnosticFileID);
end
