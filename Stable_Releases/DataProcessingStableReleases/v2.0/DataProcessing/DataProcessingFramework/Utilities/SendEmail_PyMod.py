# SendEmail_PyMod.py
#
# This Python module allows for a generic email to be composed, and has a ready-to-send email template for DP alerts
#
#   SendGenericEmail   - send email by specifying to:, from:, cc:, subject:, message:, and smtp server,login,password
#   SendDPAlert        - ready-to-send template for DP alerts. Just specify to: and message body.
#
#
# 20140717 CHF - Created
#

import smtplib

def SendGenericEmail(from_addr,to_addr_list,cc_addr_list,subject,message,login,password,smtpserver):
    """
    This function sends a generic email

    _________________________
    Inputs:
        from_addr, to_addr_list, cc_addr_list
            Email addresses, self-explanatory. Each can be list of addresses, e.g. ['john@domain.com','wayne@university.edu']
        subject, message
            Text with subject and message
        login, password
            For the smtpserver
        smtpserver
            Including port. When in doubt, use smtp.gmail.com:587 (it will require web login to authorize 3rd party emails)

    Returns:
        None

    Versioning:
        20140717 CHF - Created

    To do:
        Add return value if successful

    _________________________
    """

    # Populate message header
    header  = 'From: %s\n' % from_addr

    if to_addr_list.__class__ is list:
        header += 'To: %s\n' % ','.join(to_addr_list)
    else:
        header += 'To: %s\n' % to_addr_list

    if cc_addr_list.__class__ is list:
        header += 'Cc: %s\n' % ','.join(cc_addr_list)
    else:
        header += 'Cc: %s\n' % cc_addr_list

    header += 'Subject: %s\n\n' % subject
    message = header + message
 
    # Server login
    server = smtplib.SMTP(smtpserver)
    server.starttls()
    server.login(login,password)

    print header
    print message

    # Send email
    problems = server.sendmail(from_addr, to_addr_list, message)
    
    # Break server connection
    server.quit()

    return problems

def SendDPAlert(to_addr_list,message=''):
    """
    Pre-made email template for DP Alerts. Using ednotifier@gmail account.

    _________________________
    Inputs:
        from_addr, to_addr_list, cc_addr_list
            Email addresses, self-explanatory. Each can be list of addresses, e.g. ['john@domain.com','wayne@university.edu']
        subject, message
            Text with subject and message
        login, password
            For the smtpserver
        smtpserver
            Including port. When in doubt, use smtp.gmail.com:587 (it will require web login to authorize 3rd party emails)

    Returns:
        None

    Versioning:
        20140717 CHF - Created

    To do:
        Add return value if successful

    _________________________
    """


    from_addr = 'ednotifier@gmail.com'
    login = 'ednotifier'
    password = 'LUXalert'

    subject = 'DP alert!'

    out = SendGenericEmail(from_addr,to_addr_list,'',subject,message,login,password,'smtp.gmail.com:587')
    return out


