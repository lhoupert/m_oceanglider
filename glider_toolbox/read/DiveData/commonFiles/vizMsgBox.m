function vizMsgBox(strMsg, strTitle, strMode)
	global gDisablePopups;

    if (gDisablePopups == 1)
        disp(strMsg);
    else
        uiwait(msgbox(strMsg, strTitle, strMode));
    end
end
