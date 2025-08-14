# Use an official Python image
FROM python:3.12-slim

# Set working directory
WORKDIR /app

# Copy API code into container
COPY . /app

# Install dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Expose API port
EXPOSE 8000

# Default command to run the API
CMD ["uvicorn", "api:app", "--host", "0.0.0.0", "--port", "8000"]
