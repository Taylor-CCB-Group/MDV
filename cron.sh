!/bin/bash

# Define variables
LOG_FILE="/home/jhire/mdv_work/update_containers.log"
ERROR_LOG_FILE="/home/jhire/mdv_work/update_containers_error.log"
CURRENT_TIME=$(date '+%Y-%m-%d %H:%M:%S')
MAIL_RECIPIENT="jayesh.hire@well.ox.ac.uk"
MAIL_SUBJECT="Test Cron Job Failed"
MAIL_BODY="The test cron job failed. Please check the attached error log."

# Ensure log files exist or create them
touch "$LOG_FILE" "$ERROR_LOG_FILE"

{
echo "=================================================================================="
echo "$CURRENT_TIME"

# Change to the directory containing your Docker Compose file
echo "Changing directory to $COMPOSE_DIR..."
if ! cd "$COMPOSE_DIR"; then
    echo "==================================================================================" >> "$ERROR_LOG_FILE"
    echo "$CURRENT_TIME" >> "$ERROR_LOG_FILE"
    echo "Simulating an error to test email sending" >> "$ERROR_LOG_FILE"
fi

} >> "$LOG_FILE" 2>> "$ERROR_LOG_FILE"

# Send mail if there are errors
if [ -s "$ERROR_LOG_FILE" ]; then
echo "$MAIL_BODY" | mail -s "$MAIL_SUBJECT" "$MAIL_RECIPIENT"
echo "Email sent to $MAIL_RECIPIENT" # Debugging line
fi