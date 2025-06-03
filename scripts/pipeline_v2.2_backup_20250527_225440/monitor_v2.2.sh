#!/bin/bash

# Monitor RNA-seq Pipeline v2.2 Progress

echo "=== RNA-SEQ PIPELINE v2.2 MONITOR ==="
echo "Time: $(date)"
echo ""

# Check screen session
echo "=== SCREEN SESSION STATUS ==="
SCREEN_STATUS=$(screen -list | grep rnaseq_v2.2 || echo "NOT RUNNING")
echo "$SCREEN_STATUS"
echo ""

# Find the most recent log file
LATEST_LOG=$(ls -t /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2/pipeline_v2.2_run_*.log 2>/dev/null | head -1)

if [ -n "$LATEST_LOG" ]; then
    echo "=== LATEST LOG FILE ==="
    echo "Log file: $LATEST_LOG"
    echo "Log size: $(ls -lh "$LATEST_LOG" | awk '{print $5}')"
    echo ""
    
    # Check if pipeline is still running
    if pgrep -f "run_rnaseq_pipeline.sh" > /dev/null; then
        echo "STATUS: PIPELINE IS RUNNING ✅"
    else
        echo "STATUS: PIPELINE COMPLETED OR STOPPED ⚠️"
    fi
    
    echo ""
    echo "=== CURRENT STAGE ==="
    # Find the most recent stage message
    CURRENT_STAGE=$(grep "Stage.*:" "$LATEST_LOG" | tail -1)
    if [ -n "$CURRENT_STAGE" ]; then
        echo "$CURRENT_STAGE"
    else
        echo "No stage information found"
    fi
    
    echo ""
    echo "=== RECENT ACTIVITY (last 15 lines) ==="
    tail -15 "$LATEST_LOG"
    
    echo ""
    echo "=== DATASET PROGRESS ==="
    # Count completed datasets
    COMPLETED_DATASETS=$(grep -i "successfully processed.*data" "$LATEST_LOG" | wc -l)
    echo "Completed dataset processing steps: $COMPLETED_DATASETS"
    
    # Show dataset completion messages
    echo "Recent completions:"
    grep -i "successfully processed.*data" "$LATEST_LOG" | tail -3
    
    echo ""
    echo "=== ERROR CHECK ==="
    # Check for recent errors
    RECENT_ERRORS=$(tail -100 "$LATEST_LOG" | grep -i "error\|failed" | wc -l)
    if [ $RECENT_ERRORS -gt 0 ]; then
        echo "⚠️  Recent errors/failures found: $RECENT_ERRORS"
        echo "Most recent errors:"
        tail -100 "$LATEST_LOG" | grep -i "error\|failed" | tail -3
    else
        echo "✅ No recent errors detected"
    fi
    
    echo ""
    echo "=== V2.2 FIX STATUS ==="
    # Check for validation messages
    if grep -q "MAGE.*tissue.*mapping\|tissue.*validation" "$LATEST_LOG" 2>/dev/null; then
        echo "✅ MAGE tissue mapping processing detected"
    else
        echo "⏳ MAGE tissue mapping not yet reached"
    fi
    
    if grep -q "ADNI.*data_type\|data_type.*validation" "$LATEST_LOG" 2>/dev/null; then
        echo "✅ ADNI data_type processing detected"
    else
        echo "⏳ ADNI data_type validation not yet reached"
    fi
    
    if grep -q "ENCODE.*gene.*format\|gene_id_format" "$LATEST_LOG" 2>/dev/null; then
        echo "✅ ENCODE gene ID format processing detected"
    else
        echo "⏳ ENCODE gene ID format validation not yet reached"
    fi
    
else
    echo "ERROR: No log files found for v2.2 pipeline"
fi

echo ""
echo "=== NEXT CHECK ==="
echo "Run: ./monitor_v2.2.sh"
echo "Or attach to screen: screen -r rnaseq_v2.2"
echo "Or check logs: tail -f $LATEST_LOG"