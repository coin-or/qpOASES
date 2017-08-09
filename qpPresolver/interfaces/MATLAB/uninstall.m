%% uninstall
%  Remove all compiled files.
function uninstall

cmdline = 'make purge';
system(cmdline);

end
