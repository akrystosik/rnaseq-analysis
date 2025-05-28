#!/bin/bash

LOG_FILE="pipeline_run_20250523_233507.log"

echo "=== RNA-SEQ PIPELINE MONITOR ==="
echo "Time: $(date)"
echo "Screen session: $(screen -list | grep rnaseq_pipeline || echo 'NOT RUNNING')"

if [ -f "$LOG_FILE" ]; then
    echo "Log file size: $(ls -lh $LOG_FILE | awk '{print $5}')"
    echo ""
    
    # Check if pipeline is still running
    if pgrep -f "run_rnaseq_pipeline.sh" > /dev/null; then
        echo "STATUS: PIPELINE IS RUNNING"
    else
        echo "STATUS: PIPELINE COMPLETED OR STOPPED"
    fi
    
    echo ""
    echo "=== CURRENT STAGE ==="
    # Find the most recent stage message
    CURRENT_STAGE=$(grep "Stage.*:" $LOG_FILE | tail -1)
    if [ -n "$CURRENT_STAGE" ]; then
        echo "$CURRENT_STAGE"
    else
        echo "No stage information found"
    fi
    
    echo ""
    echo "=== RECENT ACTIVITY (last 10 lines) ==="
    tail -10 $LOG_FILE
    
    echo ""
    echo "=== DATASET PROGRESS ==="
    # Count completed datasets
    COMPLETED_DATASETS=$(grep "Successfully processed.*data:" $LOG_FILE | wc -l)
    echo "Completed datasets: $COMPLETED_DATASETS"
    
    # Show dataset completion messages
    grep "Successfully processed.*data:" $LOG_FILE | tail -5
    
    echo ""
    echo "=== ERROR CHECK ==="
    # Check for recent errors
    RECENT_ERRORS=$(tail -100 $LOG_FILE | grep -i "error\|failed" | wc -l)
    if [ $RECENT_ERRORS -gt 0 ]; then
        echo "Recent errors/failures found: $RECENT_ERRORS"
        tail -100 $LOG_FILE | grep -i "error\|failed" | tail -3
    else
        echo "No recent errors detected"
    fi
    
else
    echo "ERROR: Log file $LOG_FILE not found"
fi

echo ""
echo "=== NEXT CHECK ==="
echo "Run: ./monitor_pipeline.sh"
echo "Or attach to screen: screen -r rnaseq_pipeline"