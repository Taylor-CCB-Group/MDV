# Define variables
COMPOSE_DIR="/home/jhire/ddd"
LOG_FILE="/home/jhire/mdv_work/update_containers.log"
ERROR_LOG_FILE="/home/jhire/mdv_work/update_containers_error.log"
CURRENT_TIME=$(date '+%Y-%m-%d %H:%M:%S')
MAIL_RECIPIENT="jayesh.hire@well.ox.ac.uk"
MAIL_SUBJECT="Test Cron Job Failed"
MAIL_BODY="The test cron job failed. Please check the attached error log."

# Ensure log files exist or create them
touch "$LOG_FILE" "$ERROR_LOG_FILE"

TEMP_ERROR_LOG=$(mktemp)

echo "*********"

echo "Simulating an error to test email sending" >> "$TEMP_ERROR_LOG"

log_error() {
  echo "==================================================================================" >> "$TEMP_ERROR_LOG"
  echo "$CURRENT_TIME££££££££££££" >> "$TEMP_ERROR_LOG"
  echo "$1" >> "$TEMP_ERROR_LOG"
  echo "Error logged: $1"
}

{
echo "=================================================================================="
echo "$CURRENT_TIME"

# Change to the directory containing your Docker Compose file
echo "Changing directory to $COMPOSE_DIR..."


echo "Simulating Docker command failure..."
  if ! docker compose -f "$COMPOSE_FILE" up -d --invalid-flag; then
    log_error "Failed: Docker command failed due to invalid flag"
  fi



echo "Update completed."

} >> "$LOG_FILE" 2>> "$TEMP_ERROR_LOG"

# Debugging messages
if [ -e "$TEMP_ERROR_LOG" ]; then
  echo "Temporary error log file exists."  # Debugging line
  FILE_SIZE=$(stat -c%s "$TEMP_ERROR_LOG")
  cat "$TEMP_ERROR_LOG"
  echo "Temporary error log file size: $FILE_SIZE bytes"  # Debugging line
fi  

# Check if the temporary error log file is empty or not
if [ -s "$TEMP_ERROR_LOG" ]; then
  echo "Error log file is not empty. Preparing to send email."
  EMAIL_BODY="The cron job failed. Please find the details below:\n\n"
  EMAIL_BODY+="Timestamp: $CURRENT_TIME\n\n"
  EMAIL_BODY+="Error Details:\n"
  EMAIL_BODY+=$(cat "$TEMP_ERROR_LOG")

  echo -e "$EMAIL_BODY" | mail -s "$MAIL_SUBJECT" "$MAIL_RECIPIENT"
  echo "Email sent to $MAIL_RECIPIENT"  # Debugging line

  # Append new errors to the main error log file
  cat "$TEMP_ERROR_LOG" >> "$ERROR_LOG_FILE"
else
  echo "Error log file is empty. No email sent."  # Debugging line
fi

# Clean up the temporary error log file
rm -f "$TEMP_ERROR_LOG"